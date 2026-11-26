#!/usr/bin/env python3
"""Step 13 - isoform-specific RT-PCR / qPCR primers for nesprin-1alpha2 in human, mouse and rat.

Nesprin-1alpha2 is transcribed from an internal SYNE1/Syne1 promoter: its unique first exon
(exon 1alpha2) splices to the first shared exon (the point where alpha2 joins the giant reading
frame). The giant isoform never splices in exon 1alpha2 (it is intronic there), so a primer pair
with the FORWARD primer in exon 1alpha2 (or spanning the exon 1alpha2 -> shared-exon junction) and
the REVERSE primer in a downstream shared exon amplifies alpha2 and only alpha2. Because the first
intron is several kb, the product is also cDNA-specific (genomic DNA gives no band / a huge band).

alpha2 transcripts: human NM_001347702.2, mouse NM_022027.3 (both experimentally defined
nesprin-1alpha2), rat XM_006227852.6 (predicted; the unambiguous alpha2 ortholog - 17 exons like
mouse, shared KASH C-terminus, and an exon-1alpha2 peptide MVVAEDL...RVKKLKETLV matching human/rat).

Tm: nearest-neighbour (SantaLucia 1998), Primer3-like salt (50 mM Na+, 1.5 mM Mg2+, 0.6 mM dNTP,
250 nM oligo). Forward primers are kept in the alpha2-specific sequence and verified ABSENT from the
giant mRNA. These are vetted starting points - confirm with Primer-BLAST against the species
transcriptome before ordering.
"""
from __future__ import annotations
import json
import sys
from pathlib import Path

from Bio.SeqUtils import MeltingTemp as mt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DATA, RESULTS, revcomp, save_json  # noqa: E402

SPECIES = {  # view -> (alpha2 transcript id, giant transcript id)
    "human_SYNE1":  ("NM_001347702.2", "NM_182961.4"),
    "mouse_Syne1":  ("NM_022027.3",    "NM_001347732.2"),
    "rat_Syne1":    ("XM_006227852.6", "NM_001419849.1"),
}
LABEL = {"human_SYNE1": "human (SYNE1)", "mouse_Syne1": "mouse (Syne1)", "rat_Syne1": "rat (Syne1)"}


def genomic(view: str) -> str:
    lines = (DATA / "genomic" / f"{view}.fa").read_text().splitlines()
    return "".join(l.strip() for l in lines[1:] if not l.startswith(">")).upper()


def build_mrna(view: str, txid: str):
    """Return (mrna 5'->3', exon_bounds[(s,e) in mRNA], cds_start_idx, first_intron_bp)."""
    info = json.load(open(DATA / "models_built.json"))[view]
    st = info["region"]["start"]
    g = genomic(view)
    t = [x for x in info["transcripts"] if x["id"] == txid][0]
    ex, strand, cds = t["exons"], t["strand"], t["cds"]
    order = ex if ex[0][0] > ex[-1][0] else list(reversed(ex))  # 5'->3' transcript order
    if strand == -1:
        mrna = "".join(revcomp(g[a - st:b - st + 1]) for a, b in order)
    else:
        mrna = "".join(g[a - st:b - st + 1] for a, b in order)
    bounds, acc = [], 0
    for a, b in order:
        bounds.append((acc, acc + (b - a + 1))); acc += (b - a + 1)
    atg_g = cds[1] if strand == -1 else cds[0]
    cds_start, pos = None, 0
    for (a, b) in order:
        if a <= atg_g <= b:
            cds_start = pos + ((b - atg_g) if strand == -1 else (atg_g - a)); break
        pos += (b - a + 1)
    e1, e2 = order[0], order[1]
    intron1 = (e1[0] - e2[1] - 1) if strand == -1 else (e2[0] - e1[1] - 1)
    return mrna, bounds, cds_start, intron1


def tm(s): return round(mt.Tm_NN(s, Na=50, Mg=1.5, dNTPs=0.6, dnac1=250, dnac2=0, saltcorr=7), 1)
def gc(s): return round(100 * (s.count("G") + s.count("C")) / len(s), 1)


def ok(s):
    if not 40 <= gc(s) <= 62:
        return False
    if any(b * 5 in s for b in "ACGT"):
        return False
    if s[-1] not in "GC" or s[-5:].count("G") + s[-5:].count("C") > 3:   # 3' clamp, not too GC-rich
        return False
    return 58.0 <= tm(s) <= 63.0


def self3(s, k=4): return revcomp(s[-k:]) in s


def design(view, a2, giant):
    mrna, b, cds_start, intron1 = build_mrna(view, a2)
    gmrna, _, _, _ = build_mrna(view, giant)              # giant mRNA, for specificity check
    e1_end = b[0][1]                                       # end of exon 1alpha2 in mRNA
    # forward primers: alpha2-specific. Start from 6 nt before the ATG (keep in alpha2 coding where
    # possible) up to the junction. Junction primers need >=12 nt in exon1a2 and 4-10 nt in exon2.
    lo = max(0, cds_start - 6)
    fwd = []
    for L in range(19, 25):
        for i in range(lo, e1_end):
            s = mrna[i:i + L]
            into2 = i + L - e1_end
            if into2 > 0 and (into2 < 4 or into2 > 10 or e1_end - i < 12):
                continue
            if not ok(s) or self3(s):
                continue
            if s in gmrna or revcomp(s) in gmrna:          # must be absent from the giant isoform
                continue
            fwd.append(dict(seq=s, start=i + 1, end=i + L, tm=tm(s), gc=gc(s), junction=into2 > 0))
    # reverse primers: in shared exons 2-4 (exons index 1,2,3)
    rev = []
    for ex_i in (1, 2, 3):
        a, e = b[ex_i]
        for L in range(19, 25):
            for i in range(a, e - L + 1):
                s = revcomp(mrna[i:i + L])
                if ok(s) and not self3(s):
                    rev.append(dict(seq=s, start=i + 1, end=i + L, tm=tm(s), gc=gc(s), exon=ex_i + 1))
    pairs = []
    for f in fwd:
        for r in rev:
            size = r["end"] - f["start"] + 1
            if not (90 <= size <= 400 and abs(f["tm"] - r["tm"]) <= 1.5):
                continue
            if revcomp(f["seq"][-5:]) in r["seq"] or revcomp(r["seq"][-5:]) in f["seq"]:
                continue
            pen = abs(f["tm"] - 60) + abs(r["tm"] - 60) + abs(f["gc"] - 50) / 10 + abs(r["gc"] - 50) / 10
            pen += 0 if f["junction"] else 0.15
            pen += 0 if size <= 200 else 0.1               # mild preference for qPCR-sized products
            pairs.append(dict(forward=f, reverse=r, product=size, penalty=round(pen, 2),
                              qpcr_suitable=size <= 200))
    pairs.sort(key=lambda p: p["penalty"])
    n_junction_fwd = len({f["seq"] for f in fwd if f["junction"]})
    chosen, used_fwd, used_key = [], set(), set()
    # pass 1: ensure the single best junction-spanning pair is included if any exists
    for p in pairs:
        if p["forward"]["junction"]:
            chosen.append(p); used_fwd.add(p["forward"]["seq"])
            used_key.add((p["forward"]["seq"], p["reverse"]["exon"], p["product"] // 40)); break
    # pass 2: fill up to 4, preferring a second distinct forward and varied product sizes
    for p in pairs:
        k = (p["forward"]["seq"], p["reverse"]["exon"], p["product"] // 40)
        if k in used_key:
            continue
        # once we have 2 pairs, require the forward to be a new one until both forwards are used
        if len(chosen) >= 2 and len(used_fwd) < 2 and p["forward"]["seq"] in used_fwd:
            continue
        chosen.append(p); used_fwd.add(p["forward"]["seq"]); used_key.add(k)
        if len(chosen) == 4:
            break
    chosen.sort(key=lambda p: (not p["forward"]["junction"], p["penalty"]))
    return dict(view=view, transcript=a2, giant=giant, exon1a2_len=e1_end, cds_start_mrna=cds_start + 1,
                first_intron_bp=intron1, n_forward_candidates=len({f["seq"] for f in fwd}),
                n_junction_forward=n_junction_fwd, pairs=chosen)


def main():
    out, lines = {}, []
    hdr = ("species\ttranscript\tpair\tforward (5'->3')\tF_Tm\tF_GC%\tF_site\treverse (5'->3')\t"
           "R_Tm\tR_GC%\tR_exon\tproduct_bp\tqPCR?")
    lines.append(hdr)
    for view, (a2, giant) in SPECIES.items():
        r = design(view, a2, giant)
        out[view] = r
        print(f"\n### {LABEL[view]}  alpha2={a2}  (exon 1alpha2 = {r['exon1a2_len']} nt; "
              f"first intron {r['first_intron_bp']} bp; "
              f"{r['n_forward_candidates']} fwd candidates, {r['n_junction_forward']} junction-spanning)")
        for k, p in enumerate(r["pairs"], 1):
            f, rv = p["forward"], p["reverse"]
            fsite = "exon1a2/exon2 junction" if f["junction"] else "exon 1a2 (alpha2-unique)"
            lines.append(f"{LABEL[view]}\t{a2}\t{k}\t{f['seq']}\t{f['tm']}\t{f['gc']}\t{fsite}\t"
                         f"{rv['seq']}\t{rv['tm']}\t{rv['gc']}\texon {rv['exon']}\t{p['product']}\t"
                         f"{'yes' if p['qpcr_suitable'] else 'no'}")
            print(f"  pair {k}: F {f['seq']} (Tm {f['tm']}, {fsite})")
            print(f"          R {rv['seq']} (Tm {rv['tm']}, exon {rv['exon']})  -> {p['product']} bp"
                  f"{'  [qPCR-sized]' if p['qpcr_suitable'] else ''}")
    save_json(RESULTS / "mammalian_alpha2_primers.json", out)
    (RESULTS / "mammalian_alpha2_primers.tsv").write_text("\n".join(lines) + "\n")
    print(f"\nwrote {RESULTS/'mammalian_alpha2_primers.tsv'}")


if __name__ == "__main__":
    main()
