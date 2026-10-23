# work out which file belongs to each of the 20 comparisons
# usage: resolve_inputs.py base_dir
import csv
import sys
import zipfile
from pathlib import Path

base = Path(sys.argv[1])
samples = list(csv.DictReader(open(base / "configs" / "wgs_samples.csv")))
extract = base / "inputs" / "cleaned_ab"
extract.mkdir(parents=True, exist_ok=True)
zipfile.ZipFile(base / "cleaned_acinetobacter_usable.zip").extractall(extract)

lookup = {
    "LA1": "19606.gb",
    "LA2": "AbATCC_col_only_P8_Ab_cleaned.fasta",
    "LA3": "AbATCC_as2o3_P8_Ab_cleaned.fasta",
    "LA5": "AbATCC_col_as2o3_P8_Ab_cleaned.fasta",
    "LA8": "AbATCC_col_as2o3_P16_Ab_cleaned.fasta",
}

for s in samples:
    code = s["sample_code"]
    if s["background"] == "Ab ATCC 19606":
        name = lookup[code]
        path = (base if name.endswith(".gb") else extract) / name
    elif s["background"] == "Ab clinical KR6151":
        path = extract / f"{code}_assembly_Ab_cleaned.fasta"
    else:
        path = base / f"{code}_assembly.fasta"
    ok = path.exists()
    s["resolved_source"] = str(path) if ok else ""
    s["resolved"] = ok
    if code == "LA1":
        s["source_kind"] = "genbank_reference"
    else:
        s["source_kind"] = "cleaned_assembly" if s["background"].startswith("Ab") else "original_assembly"

outfile = base / "inputs" / "wgs_resolved_manifest.csv"
with open(outfile, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=samples[0].keys())
    w.writeheader()
    w.writerows(samples)
print("resolved", sum(s["resolved"] for s in samples), "of", len(samples))
