# match sample codes to assembly files and do a quick size check
# usage: wgs_inputs.py --mapping samples.csv --assemblies dir --out qc.csv
import argparse
import csv
from pathlib import Path

import pandas as pd
from Bio import SeqIO


def run(mapping, adir, out):
    rows = []
    for r in csv.DictReader(open(mapping)):
        code = r["sample_code"]
        f = Path(adir) / f"{code}_assembly.fasta"
        if not f.exists():
            # LA5 got renamed at some point
            alts = {"LA5": ["LA5_assembly_Acinetobacter_baumannii_cleaned.fasta", "FILTERED_Ab_col_as2o3_P8.fasta"]}
            found = [str(Path(adir) / x) for x in alts.get(code, []) if (Path(adir) / x).exists()]
            rows.append({**r, "path": "", "present": "NO", "alternatives": ";".join(found),
                         "bp": 0, "contigs": 0, "largest_contig": 0, "qc_note": "missing"})
            continue
        lens = [len(x.seq) for x in SeqIO.parse(f, "fasta")]
        size = sum(lens)
        if size > 8_000_000:
            note = "too big, maybe mixed/contaminated"
        elif size < 3_000_000:
            note = "short"
        else:
            note = "size ok"
        rows.append({**r, "path": str(f), "present": "YES", "alternatives": "", "bp": size,
                     "contigs": len(lens), "largest_contig": max(lens) if lens else 0, "qc_note": note})
    d = pd.DataFrame(rows)
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    d.to_csv(out, index=False)
    print(len(d), "samples,", (d.present == "YES").sum(), "with assembly,", d.qc_note.str.startswith("too big").sum(), "oversized")
    return d


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--mapping", required=True)
    ap.add_argument("--assemblies", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    run(a.mapping, a.assemblies, a.out)
