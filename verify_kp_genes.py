# check the Kp target genes (lptD, mgrB, lpxD, lpxC) against each assembly
# usage: verify_kp_genes.py ref.gbff manifest.csv outdir
# only looks at these genes, it's not a genome-wide call set
import csv
import sys
from pathlib import Path

from Bio import Align, SeqIO
from Bio.Seq import Seq

gbk, manifest_file, outdir = sys.argv[1:4]
outdir = Path(outdir)
outdir.mkdir(parents=True, exist_ok=True)

ref = next(SeqIO.parse(gbk, "genbank"))
genes = ["lptd", "mgrb", "lpxd", "lpxc"]
feats = {}
for f in ref.features:
    if f.type != "CDS" or "translation" not in f.qualifiers:
        continue
    name = f.qualifiers.get("gene", [""])[0].lower()
    prod = f.qualifiers.get("product", [""])[0].lower()
    for g in genes:
        if name == g or (not name and g in prod.split()):
            feats[g] = f
print("found:", list(feats))

aligner = Align.PairwiseAligner(mode="local", match_score=3, mismatch_score=-3,
                                open_gap_score=-7, extend_gap_score=-1)


def get_seeds(cds):
    # short genes (mgrB) need shorter seeds
    if len(cds) < 250:
        return [(i, cds[i:i + 17]) for i in range(5, len(cds) - 20, 19)]
    return [(i, cds[i:i + 35]) for i in range(45, len(cds) - 55, 90)]


def best_region(seeds, contigs):
    best = None
    for cid, fwd, rev in contigs:
        for strand, s in (("+", fwd), ("-", rev)):
            starts = []
            for off, seed in seeds:
                ix = s.find(seed)
                if ix >= 0:
                    starts.append(ix - off)
            if not starts:
                continue
            starts.sort()
            mid = starts[len(starts) // 2]
            votes = sum(abs(x - mid) < 100 for x in starts)
            if best is None or votes > best[0]:
                best = (votes, cid, strand, s, mid)
    return best


cds_seq = {g: str(f.location.extract(ref.seq)).upper() for g, f in feats.items()}
manifest = list(csv.DictReader(open(manifest_file)))
rows = []
for m in manifest:
    if not m["background"].startswith("Kp"):
        continue
    contigs = []
    for c in SeqIO.parse(m["resolved_source"], "fasta"):
        contigs.append((c.id, str(c.seq).upper(), str(c.seq.reverse_complement()).upper()))

    for g, f in feats.items():
        cds = cds_seq[g]
        refaa = f.qualifiers["translation"][0]
        row = {k: m[k] for k in ("background", "sample_code", "condition", "passage")}
        row.update(gene=g, reference_protein=f.qualifiers.get("protein_id", [""])[0],
                   reference_length=len(refaa), status="NO_HOMOLOG")
        best = best_region(get_seeds(cds), contigs)
        if best is not None and best[0] >= 2:
            votes, cid, strand, s, start = best
            lo = max(0, start - 100)
            hi = min(len(s), start + len(cds) + 100)
            cand = s[lo:hi]
            aln = aligner.align(cds, cand)[0]
            rblocks, qblocks = aln.aligned
            span = sum(e - b for b, e in rblocks)
            full = span >= len(cds) * 0.95 and rblocks[0][0] < 20 and rblocks[-1][1] > len(cds) - 20

            mapping = {}
            ins, dels = [], []
            for k, ((rs, re_), (qs, qe)) in enumerate(zip(rblocks, qblocks)):
                if k > 0:
                    prev_re = rblocks[k - 1][1]
                    prev_qe = qblocks[k - 1][1]
                    if rs > prev_re:
                        dels.append(f"{prev_re + 1}-{rs}")
                    if qs > prev_qe:
                        ins.append(f"after_{prev_re}:{cand[prev_qe:qs]}")
                for i, j in zip(range(rs, re_), range(qs, qe)):
                    mapping[i] = cand[j]

            nt = [f"{cds[i]}{i + 1}{mapping[i]}" for i in sorted(mapping) if cds[i] != mapping[i]]
            # only call aa changes where the whole codon is aligned (skips indel codons)
            aa = []
            for p in range(0, len(cds) - 2, 3):
                if all(i in mapping for i in (p, p + 1, p + 2)):
                    qc = mapping[p] + mapping[p + 1] + mapping[p + 2]
                    ra = str(Seq(cds[p:p + 3]).translate(table=11))
                    qa = str(Seq(qc).translate(table=11))
                    if ra != qa:
                        aa.append(f"{ra}{p // 3 + 1}{qa}")
            row.update(status="FULL_CDS" if full else "PARTIAL_CDS", contig=cid, orientation=strand,
                       matching_seeds=votes, aligned_ref_bases=span,
                       nt_substitutions=";".join(nt), indels=";".join(ins + dels),
                       aa_substitutions=";".join(aa))
        rows.append(row)
    print("done", m["sample_code"])

cols = list(dict.fromkeys(k for r in rows for k in r))
with open(outdir / "kp_coding_sequence_audit.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols)
    w.writeheader()
    w.writerows(rows)
print(len(rows), "rows written")
