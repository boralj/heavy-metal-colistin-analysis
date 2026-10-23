# check the reported substitutions (lpxD, pmrB, pmrC, lptD, mgrB, lpxC) in every assembly
# usage: verify_reported_proteins.py manifest.csv 19606.gb 700603.gbff outdir
# no seed hit just means we didn't find it, not that the gene is gone
# TODO: indels show up as a bunch of fake substitutions after the shift, need a real aligner
import csv
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

manifest_file, ab_gb, kp_gb, outdir = sys.argv[1:5]
outdir = Path(outdir)
outdir.mkdir(parents=True, exist_ok=True)

targets = ["lpxD", "pmrB", "pmrC", "lptD", "mgrB", "lpxC"]
by_product = {
    "pmrB": ["histidine kinase pmrb"],
    "mgrB": ["membrane protein mgrb"],
    "lptD": ["lps-assembly protein lptd", "lipopolysaccharide assembly protein lptd"],
}

feats = {}
for label, path in (("Ab", ab_gb), ("Kp", kp_gb)):
    rec = next(SeqIO.parse(path, "genbank"))
    for f in rec.features:
        if f.type != "CDS" or "translation" not in f.qualifiers:
            continue
        gene = f.qualifiers.get("gene", [""])[0].lower()
        prod = f.qualifiers.get("product", [""])[0].lower()
        for t in targets:
            if t.lower() == gene or any(p in prod for p in by_product.get(t, [])):
                feats[(label, t)] = (f.qualifiers["translation"][0], f.qualifiers.get("protein_id", [""])[0])


def locate(ref, trans):
    n = len(ref)
    anchors = [(o, ref[o:o + 20]) for o in (25, n // 5, n // 3, n // 2, 2 * n // 3, 4 * n // 5) if o + 20 < n]
    best = None
    for name, strand, frame, aa in trans:
        seen = set()
        for off, seed in anchors:
            at = aa.find(seed)
            while at >= 0:
                start = at - off
                if start >= 0 and start + n <= len(aa) and start not in seen:
                    seen.add(start)
                    q = aa[start:start + n]
                    ident = sum(a == b for a, b in zip(ref, q)) / n
                    if best is None or ident > best[0]:
                        best = (ident, name, strand, frame, start, q)
                at = aa.find(seed, at + 1)
    return best


manifest = list(csv.DictReader(open(manifest_file)))
rows = []
for m in manifest:
    src = Path(m["resolved_source"])
    if not src.exists() or src.suffix == ".gb":
        continue
    bg = "Ab" if m["background"].startswith("Ab") else "Kp"
    trans = []
    for contig in SeqIO.parse(src, "fasta"):
        s = str(contig.seq).upper()
        for strand, nt in (("+", s), ("-", str(Seq(s).reverse_complement()))):
            for frame in range(3):
                usable = (len(nt) - frame) // 3 * 3
                trans.append((contig.id, strand, frame, str(Seq(nt[frame:frame + usable]).translate(table=11))))
    for gene in targets:
        if (bg, gene) not in feats:
            continue
        ref, pid = feats[(bg, gene)]
        hit = locate(ref, trans)
        row = {k: m[k] for k in ("background", "sample_code", "condition", "passage")}
        row.update(gene=gene, ref_protein_id=pid, ref_length=len(ref))
        if hit:
            ident, ctg, strand, frame, start, q = hit
            diffs = ";".join(f"{a}{i + 1}{b}" for i, (a, b) in enumerate(zip(ref, q)) if a != b)
            row.update(status="MATCH" if ident >= 0.9 else "LOW_IDENTITY", identity=round(ident, 5),
                       contig=ctg, strand=strand, frame=frame,
                       reference_residue_117=ref[116] if len(ref) > 116 else "",
                       observed_residue_117=q[116] if len(q) > 116 else "",
                       amino_acid_changes=diffs[:12000])
        else:
            row["status"] = "UNRESOLVED_NO_SEED"
        rows.append(row)
    print("done", m["sample_code"])

cols = list(dict.fromkeys(k for r in rows for k in r))
with open(outdir / "reported_protein_sequence_audit.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, cols)
    w.writeheader()
    w.writerows(rows)
print(len(rows), "records")
