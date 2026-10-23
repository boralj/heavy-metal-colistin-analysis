# quick check of the outputs before they go out
# usage: check_release.py project_dir
import sys
from pathlib import Path

import pandas as pd

root = Path(sys.argv[1])

man = pd.read_csv(root / "inputs" / "wgs_resolved_manifest.csv")
assert len(man) == 20, len(man)
assert man.resolved.astype(str).str.lower().eq("true").all(), "unresolved assemblies"

q = pd.read_csv(root / "results" / "qpcr" / "metal_minus_colistin_all.csv")
assert len(q) == 890 and q.valid_control.sum() == 752, (len(q), q.valid_control.sum())

led = pd.read_csv(root / "results" / "wgs" / "manuscript_variant_evidence_ledger.csv")
assert len(led) == 7
assert (led.audit_status == "SUPPORTED_IN_ASSEMBLY").sum() == 3, led.audit_status.value_counts()

kp = pd.read_csv(root / "results" / "wgs" / "kp_coding_sequence_audit.csv")
assert len(kp) == 40, len(kp)
lb12 = kp[(kp.sample_code == "LB12") & (kp.gene == "lptd")]
assert lb12.aa_substitutions.fillna("").str.contains("T456P").all(), "LB12 lptD T456P missing"

print("all good")
