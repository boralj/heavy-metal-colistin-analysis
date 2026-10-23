# summarize the author log2FC values, match to colistin controls, PCA per species
# usage: qpcr_summaries.py workbook.xlsx --out outdir
# no p-values on purpose, the sheet has no replicate IDs
import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

KEY = ["Sample", "Species", "Iso_base", "Gene", "Phase", "Day"]


def colistin_only(desc):
    if not isinstance(desc, str):
        return False
    s = desc.lower()
    return "colistin" in s and not re.search(r"\+\s*(hg|cu|ag|ars)", s)


def main(workbook, out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    df = pd.read_excel(workbook, sheet_name="AllValues")
    df["Sample"] = df.Sample.astype(str).str.strip().str.lower()
    df["log2FC"] = pd.to_numeric(df.log2FC, errors="coerce")
    if df.log2FC.isna().any():
        raise ValueError("non-numeric log2FC somewhere")

    summ = (df.groupby(KEY, dropna=False)
              .agg(mean_log2FC=("log2FC", "mean"), sd_log2FC=("log2FC", "std"), n_values=("log2FC", "size"),
                   desc=("Isolate_description", "first"), Metal_label=("Metal_label", "first"))
              .reset_index())
    summ.to_csv(out / "sample_gene_summaries.csv", index=False)

    pairs = pd.read_excel(workbook, sheet_name="Corresponding", header=5).iloc[:, :2].dropna()
    pairs.columns = ["Sample", "Control"]
    pairs = pairs.astype(str).apply(lambda s: s.str.strip().str.lower())
    pairs.to_csv(out / "control_mapping.csv", index=False)

    p = summ.merge(pairs, on="Sample")
    c = summ[["Sample", "Species", "Iso_base", "Gene", "Day", "Phase", "desc", "mean_log2FC"]].rename(
        columns={"Sample": "Control", "Species": "ctrl_species", "Iso_base": "ctrl_isolate",
                 "desc": "ctrl_desc", "mean_log2FC": "ctrl_mean"})
    con = p.merge(c, on=["Control", "Gene", "Day", "Phase"], how="left", indicator=True)
    con["valid"] = ((con["_merge"] == "both") & (con.Species == con.ctrl_species)
                    & (con.Iso_base == con.ctrl_isolate) & con.ctrl_desc.apply(colistin_only))
    con["delta_log2FC"] = np.where(con.valid, con.mean_log2FC - con.ctrl_mean, np.nan)
    con.to_csv(out / "matched_contrasts.csv", index=False)

    # PCA, one profile per sample/phase/day. just for plotting
    ids = ["Sample", "Species", "Iso_base", "Phase", "Day", "Metal_label"]
    summ["Metal_label"] = summ["Metal_label"].fillna("none")
    prof = summ.pivot_table(index=ids, columns="Gene", values="mean_log2FC").reset_index()
    genes = [g for g in prof.columns if g not in ids]
    allscores = []
    for sp, sub in prof.groupby("Species"):
        X = sub[genes]
        keep = X.columns[X.notna().all()]
        if len(sub) < 3 or len(keep) < 2:
            continue
        v = X[keep].to_numpy(float)
        sd = v.std(axis=0, ddof=1)
        v = v[:, sd > 1e-10]
        sd = sd[sd > 1e-10]
        if v.shape[1] < 2:
            continue
        v = (v - v.mean(axis=0)) / sd
        pca = PCA(n_components=2).fit(v)
        z = pca.transform(v)
        sc = sub[ids].copy()
        sc["PC1"], sc["PC2"] = z[:, 0], z[:, 1]
        sc["PC1_var"], sc["PC2_var"] = pca.explained_variance_ratio_
        allscores.append(sc)
    if allscores:
        pd.concat(allscores, ignore_index=True).to_csv(out / "pca_scores.csv", index=False)

    print(len(df), "rows,", df.Sample.nunique(), "samples,", int(con.valid.sum()), "of", len(con), "contrasts usable")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("workbook")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    main(a.workbook, a.out)
