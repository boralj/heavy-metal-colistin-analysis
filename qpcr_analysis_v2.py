# qPCR: raw Ct -> log2FC (ddCt), or start from the archived log2FC sheet
# archived sheet has no replicate IDs so tests there are only exploratory
import argparse
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from scipy.stats import ttest_1samp
from sklearn.decomposition import PCA

matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = ["Species", "Iso_base", "Phase", "Day"]
GENEKEY = BASE + ["Gene"]
TESTKEY = GENEKEY + ["Metal"]


def bh(p):
    p = np.asarray(p, float)
    out = np.full(len(p), np.nan)
    ok = np.flatnonzero(np.isfinite(p))
    if len(ok) == 0:
        return out
    order = ok[np.argsort(p[ok])]
    adj = p[order] * len(order) / np.arange(1, len(order) + 1)
    out[order] = np.clip(np.minimum.accumulate(adj[::-1])[::-1], 0, 1)
    return out


def raw_ct(csv, refs, baseline="Untreated"):
    # each Sample = one biological sample, wells are averaged first
    d = pd.read_csv(csv)
    d["Ct"] = pd.to_numeric(d["Ct"], errors="coerce")
    if d.Ct.isna().any():
        raise ValueError("non-numeric Ct (undetermined wells?), deal with them first")
    key = BASE + ["Sample", "Metal"]
    m = (d.groupby(key + ["Gene"], dropna=False).Ct.agg(["mean", "size"]).reset_index()
          .rename(columns={"mean": "Ct_mean", "size": "n_technical_wells"}))
    parts = []
    for sp, g in refs.items():
        r = m[(m.Species == sp) & (m.Gene == g)]
        parts.append(r[key + ["Ct_mean"]].rename(columns={"Ct_mean": "reference_Ct"}))
    ref = pd.concat(parts, ignore_index=True)
    m = m.merge(ref, on=key, how="left", validate="many_to_one")
    if m.reference_Ct.isna().any():
        raise ValueError("some samples have no reference gene Ct")
    m["delta_Ct"] = m.Ct_mean - m.reference_Ct
    base = (m[m.Metal == baseline].groupby(GENEKEY, dropna=False).delta_Ct.mean()
             .reset_index().rename(columns={"delta_Ct": "baseline_delta_Ct"}))
    m = m.merge(base, on=GENEKEY, how="left", validate="many_to_one")
    if m.baseline_delta_Ct.isna().any():
        raise ValueError("no untreated baseline for some isolate/day/gene")
    m["delta_delta_Ct"] = m.delta_Ct - m.baseline_delta_Ct
    m["log2FC"] = -m.delta_delta_Ct
    m["fold_change"] = 2 ** m.log2FC
    return m


def archived(xlsx, sheet="AllValues"):
    d = pd.read_excel(xlsx, sheet_name=sheet)
    d["Metal"] = d["Metal"].fillna("Colistin only")
    d["log2FC"] = pd.to_numeric(d.log2FC, errors="coerce")
    return d.dropna(subset=TESTKEY + ["log2FC"]).reset_index(drop=True)


def summarize(d, do_tests=False):
    rows = []
    for k, g in d.groupby(TESTKEY, dropna=False):
        v = g.log2FC.to_numpy(float)
        row = dict(zip(TESTKEY, k))
        row.update(n_observations=len(v), mean_log2FC=v.mean(), sd_log2FC=v.std(ddof=1) if len(v) > 1 else np.nan)
        # no p-value if all values are identical
        row["p_two_sided"] = float(ttest_1samp(v, 0).pvalue) if do_tests and len(v) > 1 and v.std() > 0 else np.nan
        rows.append(row)
    t = pd.DataFrame(rows)
    t["p_BH"] = np.nan
    if do_tests:
        for _, idx in t.groupby(["Species", "Phase"], dropna=False).groups.items():
            t.loc[idx, "p_BH"] = bh(t.loc[idx, "p_two_sided"])
    return t


def contrasts(t):
    col = (t[t.Metal == "Colistin only"][GENEKEY + ["mean_log2FC"]]
           .rename(columns={"mean_log2FC": "colistin_log2FC"}))
    o = t[t.Metal != "Colistin only"].merge(col, on=GENEKEY, how="left", validate="many_to_one")
    o["delta_log2FC"] = o.mean_log2FC - o.colistin_log2FC
    return o


def pca_plots(t, out):
    res = []
    for sp, g in t.groupby("Species"):
        idx = ["Iso_base", "Metal", "Phase", "Day"]
        w = g.pivot_table(index=idx, columns="Gene", values="mean_log2FC")
        full = w.dropna()
        if len(full) < 3 or full.shape[1] < 2:
            continue
        pca = PCA(n_components=2)
        sc = pca.fit_transform(full.to_numpy())
        f = full.reset_index()[idx].copy()
        f["PC1"], f["PC2"] = sc[:, 0], sc[:, 1]
        f["Species"] = sp
        f["PC1_var"], f["PC2_var"] = pca.explained_variance_ratio_
        f["dropped_profiles"] = len(w) - len(full)
        res.append(f)
        fig, ax = plt.subplots(figsize=(7, 5))
        for metal, z in f.groupby("Metal"):
            ax.scatter(z.PC1, z.PC2, label=metal, s=20)
        ax.legend(fontsize=7)
        ax.set(xlabel="PC1", ylabel="PC2", title=f"{sp} PCA (complete cases)")
        fig.tight_layout()
        fig.savefig(out / f"PCA_{sp}.png", dpi=130)
        plt.close(fig)
    if res:
        pd.concat(res).to_csv(out / "pca_profiles.csv", index=False)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--xlsx")
    ap.add_argument("--sheet", default="AllValues")
    ap.add_argument("--ct")
    ap.add_argument("--out", required=True)
    ap.add_argument("--tests", action="store_true", help="exploratory tests on archived data")
    ap.add_argument("--pca", action="store_true")
    a = ap.parse_args()
    if bool(a.xlsx) == bool(a.ct):
        ap.error("need exactly one of --xlsx / --ct")

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    d = archived(a.xlsx, a.sheet) if a.xlsx else raw_ct(a.ct, {"Ab": "16S rRNA", "Kp": "rpoD"})
    if a.xlsx and a.tests:
        print("careful: archived rows may not be independent replicates")
    t = summarize(d, do_tests=bool(a.ct) or a.tests)
    c = contrasts(t)
    d.to_csv(out / "observations.csv", index=False)
    t.to_csv(out / "gene_summary.csv", index=False)
    c.to_csv(out / "matched_colistin_contrasts.csv", index=False)
    missing = c[c.colistin_log2FC.isna()][TESTKEY]
    missing.to_csv(out / "missing_colistin_matches.csv", index=False)
    if a.pca:
        pca_plots(t, out)
    print(len(d), "rows,", len(t), "groups,", len(missing), "without a colistin match")
