import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from qpcr_analysis_v2 import raw_ct
from wgs_variants_annot import annotate, parse_cs


def test_ct():
    rows = []
    for sample, cond, ref, target in [("u1", "Untreated", 20, 25), ("u2", "Untreated", 21, 26),
                                      ("t1", "AgNO3", 20, 23), ("t2", "AgNO3", 21, 24)]:
        for gene, ct in [("rpoD", ref), ("target", target)]:
            for tech in (0, 0.2):
                rows.append(dict(Species="Kp", Iso_base="KR6061", Phase="Peak", Day=8, Sample=sample,
                                 Gene=gene, Metal=cond, Ct=ct + tech))
    with tempfile.TemporaryDirectory() as td:
        f = Path(td) / "ct.csv"
        pd.DataFrame(rows).to_csv(f, index=False)
        d = raw_ct(f, {"Kp": "rpoD"})
    t = d[(d.Metal == "AgNO3") & (d.Gene == "target")]
    assert np.allclose(t.delta_Ct, 3)
    assert np.allclose(t.delta_delta_Ct, -2)
    assert np.allclose(t.log2FC, 2)
    assert np.allclose(t.fold_change, 4)
    assert (t.n_technical_wells == 2).all()


def test_snp():
    with tempfile.TemporaryDirectory() as td:
        p = Path(td) / "r.gb"
        rec = SeqRecord(Seq("ATGGAATTT"), id="ref", name="ref", description="test")
        rec.annotations["molecule_type"] = "DNA"
        rec.features = [SeqFeature(FeatureLocation(0, 9, strand=1), type="CDS",
                                   qualifiers={"gene": ["abc"], "transl_table": ["11"]})]
        SeqIO.write([rec], p, "genbank")
        r = [{"reference": "ref", "position_1based": 4, "ref": "g", "alt": "a", "type": "SNP"}]
        annotate(r, p)
    assert r[0]["aa_change"] == "E2K" and r[0]["consequence"] == "missense", r


def test_cs():
    v = parse_cs(":3*ag:2+tt:1-cg:4", 10, 0, "+", "ref", "query")
    assert [(x["type"], x["position_1based"]) for x in v] == [("SNP", 14), ("insertion", 16), ("deletion", 18)]


if __name__ == "__main__":
    test_ct()
    test_snp()
    test_cs()
    print("ok")
