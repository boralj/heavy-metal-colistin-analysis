# check which qPCR samples are paired with a proper colistin-only control
# usage: qpcr_control_check.py workbook.xlsx outdir
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

book, out = sys.argv[1], Path(sys.argv[2])
out.mkdir(parents=True, exist_ok=True)

d = pd.read_excel(book, sheet_name="AllValues")
pairs = pd.read_excel(book, sheet_name="Corresponding", header=5).iloc[:, :2].dropna()
pairs.columns = ["metal_sample", "control_sample"]
for c in pairs:
    pairs[c] = pairs[c].astype(str).str.lower().str.strip()
d["Sample"] = d["Sample"].astype(str).str.lower().str.strip()


def colistin_only(desc):
    # ExposureType can't be trusted (colistin-only samples are marked Metal), so use the description
    if not isinstance(desc, str):
        return False
    d_ = desc.lower()
    return "colistin" in d_ and not re.search(r"\+\s*(hg|cu|ag|ars)", d_)


meta = d.groupby("Sample").agg(description=("Isolate_description", "first")).reset_index()
meta["colistin_only"] = meta.description.apply(colistin_only)

chk = pairs.merge(meta.rename(columns={"Sample": "control_sample", "description": "control_description",
                                       "colistin_only": "control_ok"}), on="control_sample", how="left")
chk.to_csv(out / "sample_pair_validation.csv", index=False)

key = ["Sample", "Gene", "Phase", "Day"]
a = (d.groupby(key, dropna=False)
      .agg(mean_log2FC=("log2FC", "mean"), n=("log2FC", "count"), Species=("Species", "first"),
           Iso_base=("Iso_base", "first"), description=("Isolate_description", "first"))
      .reset_index())
j = a.merge(pairs, left_on="Sample", right_on="metal_sample")
ctrl = a.rename(columns={"Sample": "control_sample", "mean_log2FC": "ctrl_mean", "n": "ctrl_n",
                         "Species": "ctrl_species", "Iso_base": "ctrl_isolate", "description": "ctrl_desc"})
z = j.merge(ctrl[["control_sample", "Gene", "Phase", "Day", "ctrl_mean", "ctrl_n", "ctrl_species",
                  "ctrl_isolate", "ctrl_desc"]],
            on=["control_sample", "Gene", "Phase", "Day"], how="left", indicator=True)

z["same_isolate"] = (z.Species == z.ctrl_species) & (z.Iso_base == z.ctrl_isolate)
z["ctrl_ok"] = z.ctrl_desc.apply(colistin_only)
z["usable"] = (z["_merge"] == "both") & z.same_isolate & z.ctrl_ok
z["delta_log2FC"] = np.where(z.usable, z.mean_log2FC - z.ctrl_mean, np.nan)
z.to_csv(out / "comparison_audit.csv", index=False)

print("total", len(z), "| usable", z.usable.sum(), "| no control for gene/phase/day",
      (z["_merge"] == "left_only").sum(), "| bad control label", (z.same_isolate & ~z.ctrl_ok).sum())
print("controls that are colistin-only:", meta[meta.colistin_only & meta.Sample.isin(pairs.control_sample)].shape[0])
