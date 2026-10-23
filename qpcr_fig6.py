# fig 6 stuff: gene means, metal minus colistin deltas, PCA, heatmaps
# usage: qpcr_fig6.py workbook.xlsx outdir
# workbook only has log2FC (no replicate ids) so no p-values / volcano
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

book, out = sys.argv[1], Path(sys.argv[2])
out.mkdir(parents=True, exist_ok=True)

vals = pd.read_excel(book, sheet_name="AllValues")
pairs = pd.read_excel(book, sheet_name="Corresponding", header=5).iloc[:, :2].dropna()
pairs.columns = ["Sample", "Control"]
pairs = pairs.astype(str).apply(lambda s: s.str.strip().str.lower()).drop_duplicates()

vals["Sample"] = vals.Sample.astype(str).str.strip().str.lower()
vals["log2FC"] = pd.to_numeric(vals.log2FC, errors="coerce")
vals = vals.dropna(subset=["log2FC", "Gene", "Sample", "Phase", "Day"])
vals["Metal_label"] = vals.Metal_label.fillna("none")  # colistin-only has no label, pivot drops NaN otherwise

key = ["Sample", "Gene", "Phase", "Day"]
st = (vals.groupby(key)
          .agg(log2FC=("log2FC", "mean"), sd=("log2FC", "std"), n=("log2FC", "size"),
               Species=("Species", "first"), Isolate=("Iso_base", "first"),
               Metal=("Metal_label", "first"), desc=("Isolate_description", "first"))
          .reset_index())
st.to_csv(out / "gene_level_means.csv", index=False)

# pair each metal sample with its colistin control
pr = st.merge(pairs, on="Sample")
ctrl = st.rename(columns={"Sample": "Control", "log2FC": "ctrl_log2FC", "Species": "ctrl_species",
                          "Isolate": "ctrl_isolate", "desc": "ctrl_desc"})
pr = pr.merge(ctrl[["Control", "Gene", "Phase", "Day", "ctrl_log2FC", "ctrl_species", "ctrl_isolate", "ctrl_desc"]],
              on=["Control", "Gene", "Phase", "Day"], how="left")

d = pr.ctrl_desc.fillna("").str.lower()
col_only = d.str.contains("colistin") & ~d.str.contains(r"\+\s*(?:hg|cu|ag|ars)")
pr["valid_control"] = (pr.ctrl_log2FC.notna() & (pr.Species == pr.ctrl_species)
                       & (pr.Isolate == pr.ctrl_isolate) & col_only)
pr["delta_log2FC"] = np.where(pr.valid_control, pr.log2FC - pr.ctrl_log2FC, np.nan)
pr.to_csv(out / "metal_minus_colistin_all.csv", index=False)
pr[~pr.valid_control].to_csv(out / "unmatched_or_invalid_controls.csv", index=False)

# PCA per species
meta = ["Sample", "Species", "Isolate", "Metal", "Phase", "Day"]
prof = st.pivot_table(index=meta, columns="Gene", values="log2FC")
scores = []
for sp, one in prof.groupby(level="Species"):
    x = one.loc[:, one.notna().mean() >= 0.8].dropna()
    x = x.loc[:, x.std() > 1e-12]
    if len(x) < 3 or x.shape[1] < 2:
        continue
    m = x.to_numpy(float)
    m = (m - m.mean(axis=0)) / m.std(axis=0, ddof=1)
    pca = PCA(n_components=2)
    xy = pca.fit_transform(m)
    f = x.index.to_frame(index=False)
    f["PC1"], f["PC2"] = xy[:, 0], xy[:, 1]
    f["PC1_var"], f["PC2_var"] = pca.explained_variance_ratio_
    f["n_genes"] = x.shape[1]
    f["dropped_profiles"] = len(one) - len(x)
    scores.append(f)

    fig, ax = plt.subplots(figsize=(8, 6))
    for cond, part in f.groupby("Metal"):
        ax.scatter(part.PC1, part.PC2, label=cond, s=24, alpha=0.85)
    ax.set(xlabel=f"PC1 ({pca.explained_variance_ratio_[0] * 100:.1f}%)",
           ylabel=f"PC2 ({pca.explained_variance_ratio_[1] * 100:.1f}%)", title=sp)
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(out / f"PCA_{sp}.png", dpi=200)
    plt.close(fig)
if scores:
    pd.concat(scores, ignore_index=True).to_csv(out / "PCA_scores.csv", index=False)

# heatmaps of metal-minus-colistin
good = pr[pr.valid_control]
for (sp, iso), g in good.groupby(["Species", "Isolate"]):
    mat = g.pivot_table(index="Gene", columns=["Metal", "Phase"], values="delta_log2FC")
    if mat.empty:
        continue
    mat.to_csv(out / f"heatmap_values_{sp}_{iso}.csv")
    v = mat.to_numpy(float)
    lim = np.nanmax(np.abs(v)) if np.isfinite(v).any() else 1
    fig, ax = plt.subplots(figsize=(max(7, mat.shape[1] * 0.45), max(4, mat.shape[0] * 0.23)))
    im = ax.imshow(v, aspect="auto", cmap="coolwarm", vmin=-lim, vmax=lim)
    ax.set_yticks(range(len(mat.index)), labels=mat.index, fontsize=6)
    ax.set_xticks(range(len(mat.columns)), labels=[" ".join(map(str, k)) for k in mat.columns],
                  rotation=90, fontsize=6)
    fig.colorbar(im, ax=ax, label="metal+colistin minus colistin (log2FC)")
    fig.tight_layout()
    fig.savefig(out / f"heatmap_{sp}_{iso}.png", dpi=180)
    plt.close(fig)

print(len(vals), "obs,", len(st), "means,", int(pr.valid_control.sum()), "/", len(pr), "valid comparisons,",
      len(scores), "PCAs")
