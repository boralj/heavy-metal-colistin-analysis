# line up the variants reported in the manuscript with what the assemblies show
# usage: variant_ledger.py results_dir
import sys
from pathlib import Path

import pandas as pd

res = Path(sys.argv[1])
ab = pd.read_csv(res / "ab_atcc_protein_audit.csv").fillna("")
kp = pd.read_csv(res / "kp_coding_sequence_audit.csv").fillna("")

claims = [
    ("Ab ATCC 19606", "LA5", "lpxD", "E117K", "reported at P8"),
    ("Ab ATCC 19606", "LA8", "lpxD", "E117K", "reported retained at P16"),
    ("Ab ATCC 19606", "LA5", "pmrB", "N440H", "reported at P8"),
    ("Ab ATCC 19606", "LA8", "pmrB", "N440H", "reported retained at P16"),
    ("Ab ATCC 19606", "LA5", "pmrB", "A444F", "reported at P8"),
    ("Ab ATCC 19606", "LA8", "pmrB", "A444F", "reported retained at P16"),
    ("Kp ATCC 700603", "LB12", "lptD", "T456P", "reported P16"),
]

out = []
for bg, sample, gene, change, claim in claims:
    if bg.startswith("Ab"):
        hits = ab[(ab.sample_code == sample) & ab.gene.str.contains(gene.lower(), case=False, regex=False)]
        col = "changes"
    else:
        hits = kp[(kp.sample_code == sample) & (kp.gene.str.lower() == gene.lower())]
        col = "aa_substitutions"

    if hits.empty:
        status, details = "NOT_EVALUATED", "no alignment record"
    else:
        good = hits[hits.status.isin(["ALIGNED", "FULL_CDS"])]
        if good.empty:
            status, details = "UNRESOLVED", "alignment incomplete or missing"
        else:
            found = good[col].astype(str).apply(lambda s: change in s.split(";"))
            status = "SUPPORTED_IN_ASSEMBLY" if found.any() else "NOT_OBSERVED_IN_ASSEMBLY"
            details = " | ".join(good[col].astype(str)) or "none"
    out.append(dict(background=bg, sample_code=sample, gene=gene, reported_change=change,
                    manuscript_claim=claim, audit_status=status, observed_changes=details))

df = pd.DataFrame(out)
df.to_csv(res / "manuscript_variant_evidence_ledger.csv", index=False)
print(df[["sample_code", "gene", "reported_change", "audit_status"]].to_string(index=False))
