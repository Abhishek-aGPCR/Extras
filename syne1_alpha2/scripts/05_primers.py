#!/usr/bin/env python3
"""Step 5 - in-silico RT-PCR primer suggestions for validating the predicted zebrafish
syne1a alpha2 transcript (GRCz11 model, +DV23).

Forward primers lie in the unique exon 1alpha2 (or span its junction with the first shared
exon); reverse primers lie in the 2nd-4th exons of the alpha2 model. Products therefore span
introns (cDNA-specific) and cannot come from the giant isoform. Tm: nearest-neighbour
(SantaLucia 1998) with Primer3-like salt settings (50 mM Na+, 1.5 mM Mg2+, 0.6 mM dNTP,
250 nM oligo). These are starting points: check specificity with Primer-BLAST before ordering.
"""
from __future__ import annotations

import sys
from pathlib import Path

from Bio.SeqUtils import MeltingTemp as mt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DATA, RESULTS, load_json, revcomp, save_json  # noqa: E402


def tm(s: str) -> float:
    return round(mt.Tm_NN(s, Na=50, Mg=1.5, dNTPs=0.6, dnac1=250, dnac2=0, saltcorr=7), 1)


def gc(s: str) -> float:
    return round(100 * (s.count("G") + s.count("C")) / len(s), 1)


def ok(s: str) -> bool:
    if not 40 <= gc(s) <= 62:
        return False
    if any(b * 5 in s for b in "ACGT"):
        return False
    if s[-1] not in "GC" or s[-5:].count("G") + s[-5:].count("C") > 3:
        return False
    return 58.0 <= tm(s) <= 63.0


def self_compl_3p(s: str, k: int = 4) -> bool:
    """3' end complementary to any part of itself (primer-dimer risk)."""
    return revcomp(s[-k:]) in s


def main() -> None:
    built = load_json(DATA / "models_built.json")["zebrafish_syne1a"]
    t = [x for x in built["transcripts"] if x["id"] == "syne1a-alpha2-pred-plusDV23-GRCz11"][0]
    st = built["region"]["start"]
    g = "".join(open(DATA / "genomic" / "zebrafish_syne1a.fa").read().split("\n")[1:])
    exons = [g[a - st:b - st + 1] for a, b in t["exons"]]
    bounds, acc = [], 0
    for e in exons:
        bounds.append((acc, acc + len(e)))
        acc += len(e)
    mrna = "".join(exons)
    e1 = bounds[0]
    utr = t["cds"][0] - t["exons"][0][0]                 # 5'UTR length in exon 1alpha2 (approximate TSS)
    fwd = []
    for L in range(19, 25):
        # forward primers start no more than 5 nt upstream of the ATG (the predicted 5'UTR/TSS is
        # uncertain); junction primers must have >=12 nt in exon 1alpha2 and >=5 nt in exon 2
        for i in range(max(0, utr - 5), e1[1]):
            s = mrna[i:i + L]
            into_ex2 = i + L - e1[1]
            if into_ex2 > 0 and (into_ex2 < 5 or into_ex2 > 10 or e1[1] - i < 12):
                continue
            if not ok(s) or self_compl_3p(s):
                continue
            fwd.append(dict(seq=s, start=i + 1, end=i + L, tm=tm(s), gc=gc(s), junction=into_ex2 > 0))
    rev = []
    for ex_i in (1, 2, 3):                                 # exons 2-4 of the alpha2 model
        a, b = bounds[ex_i]
        for L in range(19, 25):
            for i in range(a, b - L + 1):
                s = revcomp(mrna[i:i + L])
                if ok(s) and not self_compl_3p(s):
                    rev.append(dict(seq=s, start=i + 1, end=i + L, tm=tm(s), gc=gc(s), exon=ex_i + 1))
    pairs = []
    for f in fwd:
        for r in rev:
            size = r["end"] - f["start"] + 1
            if 150 <= size <= 450 and abs(f["tm"] - r["tm"]) <= 1.5:
                # avoid 3' complementarity between the two primers
                if revcomp(f["seq"][-5:]) in r["seq"] or revcomp(r["seq"][-5:]) in f["seq"]:
                    continue
                pen = abs(f["tm"] - 60) + abs(r["tm"] - 60) + abs(f["gc"] - 50) / 10 + abs(r["gc"] - 50) / 10
                pen += 0 if f["junction"] else 0.2      # slight preference for junction-spanning
                pairs.append(dict(forward=f, reverse=r, product=size, penalty=round(pen, 2)))
    pairs.sort(key=lambda p: p["penalty"])
    chosen, used = [], set()
    for p in pairs:
        key = (p["forward"]["start"] // 8, p["reverse"]["exon"])
        if key in used:
            continue
        used.add(key)
        chosen.append(p)
        if len(chosen) == 4:
            break
    intron1 = t["exons"][1][0] - t["exons"][0][1] - 1
    out = dict(transcript=t["id"], note=__doc__.strip().split("\n\n")[1].replace("\n", " "),
               first_intron_bp=intron1, pairs=chosen)
    save_json(RESULTS / "zebrafish_alpha2_primer_suggestions.json", out)
    lines = ["pair\tforward (5'->3')\tF_Tm\tF_GC\tF_location\treverse (5'->3')\tR_Tm\tR_GC\tR_location\tcDNA_product_bp"]
    for k, p in enumerate(chosen, 1):
        f, r = p["forward"], p["reverse"]
        floc = "exon1a2/exon2 junction" if f["junction"] else "exon 1a2 (unique)"
        lines.append(f"{k}\t{f['seq']}\t{f['tm']}\t{f['gc']}\t{floc} (mRNA {f['start']}-{f['end']})\t{r['seq']}\t"
                     f"{r['tm']}\t{r['gc']}\texon {r['exon']} (mRNA {r['start']}-{r['end']})\t{p['product']}")
    (RESULTS / "zebrafish_alpha2_primer_suggestions.tsv").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"(the shortest genomic product would include the {intron1} bp intron after exon 1a2)")


if __name__ == "__main__":
    main()
