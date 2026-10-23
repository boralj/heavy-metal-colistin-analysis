# metal vs colistin contrasts from the qPCR workbook, using the "Corresponding" sheet
import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd


def is_col_only(x):
    if not isinstance(x, str):
        return False
    x = x.lower()
    return "colistin" in x and not re.search(r"\+\s*(hg|cu|ag|ars)", x)


def run(path, out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    df = pd.read_excel(path, sheet_name="AllValues")
    pair = pd.read_excel(path, sheet_name="Corresponding", header=5).iloc[:, :2].dropna()
    pair.columns = ["metal_sample", "control_sample"]
    for c in pair:
        pair[c] = pair[c].astype(str).str.strip().str.lower()
    df["Sample"] = df["Sample"].astype(str).str.lower().str.strip()

    key = ["Sample", "Gene", "Phase", "Day"]
    st = (df.groupby(key, dropna=False)
            .agg(mean_log2FC=("log2FC", "mean"), sd_log2FC=("log2FC", "std"), n_values=("log2FC", "count"),
                 Species=("Species", "first"), Iso_base=("Iso_base", "first"), Metal=("Metal", "first"),
                 desc=("Isolate_description", "first"))
            .reset_index())
    m = st.merge(pair, left_on="Sample", right_on="metal_sample")
    ctrl = st[key + ["mean_log2FC", "n_values", "Species", "Iso_base", "desc"]].rename(
        columns={"Sample": "control_sample", "mean_log2FC": "ctrl_mean", "n_values": "ctrl_n",
                 "Species": "ctrl_species", "Iso_base": "ctrl_isolate", "desc": "ctrl_desc"})
    z = m.merge(ctrl, on=["control_sample", "Gene", "Phase", "Day"], how="left", indicator=True)

    z["same_isolate"] = (z.Species == z.ctrl_species) & (z.Iso_base == z.ctrl_isolate)
    z["ctrl_is_colistin"] = z.ctrl_desc.apply(is_col_only)
    z["usable"] = (z["_merge"] == "both") & z.same_isolate & z.ctrl_is_colistin
    z["delta_log2FC"] = np.where(z.usable, z.mean_log2FC - z.ctrl_mean, np.nan)

    z.to_csv(out / "mapped_metal_colistin_contrasts.csv", index=False)
    z[~z.usable].to_csv(out / "unresolved_comparisons.csv", index=False)
    st.to_csv(out / "sample_gene_phase_stats.csv", index=False)
    pair.to_csv(out / "sample_to_colistin_mapping.csv", index=False)
    print(len(df), "measurements,", df.Sample.nunique(), "samples,", int(z.usable.sum()), "usable /", len(z))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("workbook")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    run(a.workbook, a.out)
