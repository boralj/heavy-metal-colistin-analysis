# assembly vs reference with minimap2, then annotate SNPs with aa changes
import argparse
import csv
import re
import subprocess
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq


def parse_cs(cs, rstart, qstart, strand, rname, qname):
    toks = re.findall(r"[:=*+\-~][A-Za-z0-9]+", cs)
    r, q = rstart, qstart
    out = []

    def add(pos, ref, alt, typ):
        out.append(dict(reference=rname, position_1based=pos, ref=ref, alt=alt, type=typ,
                        query=qname, strand=strand, query_alignment_offset=q))

    for t in toks:
        op, body = t[0], t[1:].lower()
        if op == ":":
            r += int(body)
            q += int(body)
        elif op == "=":
            r += len(body)
            q += len(body)
        elif op == "*":
            add(r + 1, body[0], body[1], "SNP")
            r += 1
            q += 1
        elif op == "+":
            add(r, "", body, "insertion")
            q += len(body)
        elif op == "-":
            add(r + 1, body, "", "deletion")
            r += len(body)
    return out


def read_paf(path, min_id=0.95, min_len=1000):
    aln, var = [], []
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        q, qs, strand, r, rs, rend = f[0], int(f[2]), f[4], f[5], int(f[7]), int(f[8])
        match, block, mapq = int(f[9]), int(f[10]), int(f[11])
        tags = {x[:2]: x[5:] for x in f[12:]}
        if tags.get("tp") == "S":
            continue
        ident = 1 - float(tags["de"]) if "de" in tags else match / block
        ok = ident >= min_id and block >= min_len
        aln.append(dict(query=q, reference=r, strand=strand, reference_start=rs, reference_end=rend,
                        aligned_block=block, identity_estimate=ident, mapping_quality=mapq, accepted=ok))
        if ok:
            var += parse_cs(tags["cs"], rs, qs, strand, r, q)
    return aln, var


def annotate(rows, gbk):
    recs = {r.id: r for r in SeqIO.parse(gbk, "genbank")}
    for row in rows:
        row.update(genes="", feature="intergenic", aa_change="", consequence="", reference_validated="")
        rec = recs.get(row["reference"])
        if rec is None:
            row["consequence"] = "reference_id_not_found"
            continue
        pos = int(row["position_1based"]) - 1
        is_ins = row["type"] == "insertion"
        hits = []
        for ft in rec.features:
            if ft.type != "CDS":
                continue
            if not (ft.location.start <= pos < ft.location.end):
                continue
            positions = list(ft.location)  # goes in transcript order, so reversed for - strand
            if pos not in positions and not (is_ins and pos - 1 in positions):
                continue
            hits.append(ft.qualifiers.get("gene", ft.qualifiers.get("locus_tag", ["unnamed"]))[0])
            if row["type"] == "SNP":
                if str(rec.seq[pos]).upper() != row["ref"].upper():
                    row["reference_validated"] = "REF_MISMATCH"
                    row["consequence"] = "reference_allele_mismatch"
                    continue
                row["reference_validated"] = "PASS"
                old = str(ft.extract(rec.seq)).upper()
                if len(old) % 3:
                    row["consequence"] = "incomplete_cds"
                    continue
                idx = positions.index(pos)
                alt = row["alt"].upper()
                if ft.location.strand == -1:
                    alt = str(Seq(alt).reverse_complement())
                new = old[:idx] + alt + old[idx + 1:]
                cod = idx // 3
                table = int(ft.qualifiers.get("transl_table", ["11"])[0])
                a1 = str(Seq(old[cod * 3:cod * 3 + 3]).translate(table=table))
                a2 = str(Seq(new[cod * 3:cod * 3 + 3]).translate(table=table))
                row["aa_change"] = f"{a1}{cod + 1}{a2}"
                if a1 == a2:
                    row["consequence"] = "synonymous"
                elif a2 == "*":
                    row["consequence"] = "stop_gained"
                elif a1 == "*":
                    row["consequence"] = "stop_lost"
                else:
                    row["consequence"] = "missense"
            else:
                size = len(row["alt"]) if is_ins else len(row["ref"])
                row["consequence"] = "possible_frameshift" if size % 3 else "inframe_indel"
        if hits:
            row["genes"] = ";".join(dict.fromkeys(hits))
            row["feature"] = "CDS"
        elif row["type"] == "SNP":
            row["reference_validated"] = "PASS" if str(rec.seq[pos]).upper() == row["ref"].upper() else "REF_MISMATCH"
    return rows


def write_csv(path, rows, fields):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref")
    ap.add_argument("--query")
    ap.add_argument("--gbk")
    ap.add_argument("--paf")
    ap.add_argument("--out", required=True)
    ap.add_argument("--min_id", type=float, default=0.95)
    ap.add_argument("--min_len", type=int, default=1000)
    a = ap.parse_args()

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    paf = a.paf
    if not paf:
        paf = out / "alignments.paf"
        with open(paf, "w") as fh:
            subprocess.run(["minimap2", "-cx", "asm5", "--cs=long", "--secondary=no", a.ref, a.query],
                           stdout=fh, check=True)
    aln, var = read_paf(paf, a.min_id, a.min_len)
    if a.gbk:
        annotate(var, a.gbk)
    write_csv(out / "alignments.csv", aln, list(aln[0]) if aln else [])
    write_csv(out / "variants.csv", var, ["reference", "position_1based", "ref", "alt", "type", "query", "strand",
                                          "query_alignment_offset", "genes", "feature", "aa_change",
                                          "consequence", "reference_validated"])
    print(len(aln), "alignments,", len(var), "variants")
