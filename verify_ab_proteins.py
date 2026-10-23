# look for lpxD / pmrB changes in the Ab ATCC 19606 assemblies
# usage: verify_ab_proteins.py 19606.gb manifest.csv outdir
# the genbank file only has protein translations so we search six-frame translations
import csv
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

gbk, manifest_file, outdir = sys.argv[1:4]
outdir = Path(outdir)
outdir.mkdir(parents=True, exist_ok=True)

rec = next(SeqIO.parse(gbk, "genbank"))
refs = {}
for f in rec.features:
    if f.type != "CDS" or "translation" not in f.qualifiers:
        continue
    prod = f.qualifiers.get("product", [""])[0].lower()
    gene = f.qualifiers.get("gene", [""])[0].lower()
    if gene == "lpxd" or "histidine kinase pmrb" in prod or "phosphoethanolamine--lipid a transferase" in prod:
        refs[f.qualifiers.get("protein_id", [""])[0]] = (f.qualifiers["translation"][0], gene, prod)

manifest = list(csv.DictReader(open(manifest_file)))
rows = []
for m in manifest:
    if m["background"] != "Ab ATCC 19606" or m["condition"] == "untreated":
        continue
    best = {pid: None for pid in refs}
    for contig in SeqIO.parse(m["resolved_source"], "fasta"):
        seq = str(contig.seq).upper()
        for strand, nt in (("+", seq), ("-", str(Seq(seq).reverse_complement()))):
            for frame in range(3):
                usable = (len(nt) - frame) // 3 * 3
                aa = str(Seq(nt[frame:frame + usable]).translate(table=11))
                for pid, (ref, gene, prod) in refs.items():
                    n = len(ref)
                    for pos in (20, n // 4, n // 2, 3 * n // 4):
                        if pos + 18 > n:
                            continue
                        ix = aa.find(ref[pos:pos + 18])
                        if ix < pos or ix - pos + n > len(aa):
                            continue
                        cand = aa[ix - pos:ix - pos + n]
                        ident = sum(a == b for a, b in zip(cand, ref)) / n
                        if best[pid] is None or ident > best[pid][0]:
                            best[pid] = (ident, cand, contig.id, strand, frame)
    for pid, (ref, gene, prod) in refs.items():
        row = {k: m[k] for k in ("sample_code", "condition", "passage")}
        row.update(protein_id=pid, gene=gene or prod, reference_length=len(ref))
        if best[pid] is None:
            row["status"] = "UNRESOLVED"
        else:
            ident, q, ctg, strand, frame = best[pid]
            diffs = [f"{a}{i + 1}{b}" for i, (a, b) in enumerate(zip(ref, q)) if a != b]
            row.update(status="ALIGNED" if ident > 0.9 else "LOW_IDENTITY", identity=ident,
                       contig=ctg, strand=strand, frame=frame, changes=";".join(diffs),
                       residue440_reference=ref[439] if len(ref) > 439 else "",
                       residue440_query=q[439] if len(q) > 439 else "")
        rows.append(row)
    print("done", m["sample_code"])

cols = list(dict.fromkeys(k for r in rows for k in r))
with open(outdir / "ab_atcc_protein_audit.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, cols)
    w.writeheader()
    w.writerows(rows)
print(len(rows), "records")
