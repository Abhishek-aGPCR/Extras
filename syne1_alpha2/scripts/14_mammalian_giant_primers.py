#!/usr/bin/env python3
"""Step 14 - giant-isoform-specific RT-PCR / qPCR primers for nesprin-1 giant in human, mouse, rat.

The matched counterpart to step 13 (alpha2 primers). alpha2 and the giant isoform converge on a
single "shared exon" (giant exon ~130); alpha2 reaches it from its unique exon 1alpha2, the giant
reaches it from the preceding giant exon ("exon 127", a conserved 167 bp exon that alpha2 skips). So
a FORWARD primer in giant exon 127 (or spanning the giant exon127 -> shared-exon junction) plus a
REVERSE primer in the shared exon / just downstream amplifies the giant and NOT alpha2 - the mirror
image of the alpha2 assay, using the same shared-exon region so the two assays can share reverse
primers and run side by side (giant vs alpha2).

giant transcripts: human NM_182961.4, mouse NM_001347732.2, rat NM_001419849.1 (all 146 exons).
Forward primers are verified ABSENT from the alpha2 mRNA; products span a multi-kb intron
(cDNA-specific) and are validated by in-silico PCR to give one giant band and no alpha2 band.
Confirm with Primer-BLAST before ordering. Tm: nearest-neighbour (SantaLucia 1998).
"""
from __future__ import annotations
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DATA, RESULTS, revcomp, save_json  # noqa: E402
from importlib import import_module
m13 = import_module("13_mammalian_alpha2_primers")
tm, gc, ok, self3, build_mrna = m13.tm, m13.gc, m13.ok, m13.self3, m13.build_mrna

SPECIES = {  # view -> (giant id, alpha2 id)
    "human_SYNE1": ("NM_182961.4",    "NM_001347702.2"),
    "mouse_Syne1": ("NM_001347732.2", "NM_022027.3"),
    "rat_Syne1":   ("NM_001419849.1", "XM_006227852.6"),
}
LABEL = {"human_SYNE1": "human (SYNE1)", "mouse_Syne1": "mouse (Syne1)", "rat_Syne1": "rat (Syne1)"}


def shared_index(view, giant, a2):
    """Index (in giant transcript order) of the first exon shared with alpha2 (the 'shared exon')."""
    info = json.load(open(DATA / "models_built.json"))[view]
    ta = [x for x in info["transcripts"] if x["id"] == a2][0]
    tg = [x for x in info["transcripts"] if x["id"] == giant][0]
    a2o = ta["exons"] if ta["exons"][0][0] > ta["exons"][-1][0] else list(reversed(ta["exons"]))
    go = tg["exons"] if tg["exons"][0][0] > tg["exons"][-1][0] else list(reversed(tg["exons"]))
    shared = tuple(a2o[1])
    S = [i for i, e in enumerate(go) if tuple(e) == shared][0]
    # genomic first intron between giant exon S-1 and shared exon S
    strand = tg["strand"]
    e_up, e_sh = go[S - 1], go[S]
    intron = (e_up[0] - e_sh[1] - 1) if strand == -1 else (e_sh[0] - e_up[1] - 1)
    return S, intron


def design(view, giant, a2):
    gm, gb, _, _ = build_mrna(view, giant)       # full giant mRNA + exon bounds
    am, _, _, _ = build_mrna(view, a2)           # alpha2 mRNA (specificity check)
    S, intron_up = shared_index(view, giant, a2)
    up_s, up_e = gb[S - 1]                        # giant exon127 (giant-specific) span in mRNA
    sh_s, sh_e = gb[S]                            # shared exon span in mRNA
    nxt_s, nxt_e = gb[S + 1]                      # next shared exon (downstream)
    # forward primers: in giant exon127 (giant-specific) or spanning exon127->shared junction
    fwd = []
    for L in range(19, 25):
        for i in range(up_s, up_e):
            s = gm[i:i + L]
            into = i + L - up_e                   # nt into the shared exon
            if into > 0 and (into < 4 or into > 10 or up_e - i < 12):
                continue
            if not ok(s) or self3(s):
                continue
            if s in am or revcomp(s) in am:        # must be absent from alpha2
                continue
            fwd.append(dict(seq=s, start=i + 1, end=i + L, tm=tm(s), gc=gc(s), junction=into > 0))
    # reverse primers: shared exon and next shared exon (can coincide with the alpha2-assay reverses)
    rev = []
    for (a, e, tag) in [(sh_s, sh_e, "shared exon"), (nxt_s, nxt_e, "shared exon +1")]:
        for L in range(19, 25):
            for i in range(a, e - L + 1):
                s = revcomp(gm[i:i + L])
                if ok(s) and not self3(s):
                    rev.append(dict(seq=s, start=i + 1, end=i + L, tm=tm(s), gc=gc(s), loc=tag))
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
            pen += 0 if size <= 200 else 0.1
            pairs.append(dict(forward=f, reverse=r, product=size, penalty=round(pen, 2),
                              qpcr_suitable=size <= 200))
    pairs.sort(key=lambda p: p["penalty"])
    chosen, used_fwd, used_key = [], set(), set()
    for p in pairs:                               # ensure one junction-spanning pair first
        if p["forward"]["junction"]:
            chosen.append(p); used_fwd.add(p["forward"]["seq"])
            used_key.add((p["forward"]["seq"], p["reverse"]["loc"], p["product"] // 40)); break
    for p in pairs:
        k = (p["forward"]["seq"], p["reverse"]["loc"], p["product"] // 40)
        if k in used_key:
            continue
        if len(chosen) >= 2 and len(used_fwd) < 2 and p["forward"]["seq"] in used_fwd:
            continue
        chosen.append(p); used_fwd.add(p["forward"]["seq"]); used_key.add(k)
        if len(chosen) == 4:
            break
    chosen.sort(key=lambda p: (not p["forward"]["junction"], p["penalty"]))
    # in-silico PCR validation on full giant (expect [size]) and full alpha2 (expect [])
    import re
    for p in chosen:
        f, r = p["forward"]["seq"], revcomp(p["reverse"]["seq"])
        pg = [rm.end() - fm.start() for fm in re.finditer(re.escape(f), gm)
              for rm in re.finditer(re.escape(r), gm) if 0 < rm.end() - fm.start() <= 1000]
        pa = [rm.end() - fm.start() for fm in re.finditer(re.escape(f), am)
              for rm in re.finditer(re.escape(r), am) if 0 < rm.end() - fm.start() <= 1000]
        p["giant_products"], p["alpha2_products"] = pg, pa
        p["validated"] = (pg == [p["product"]] and pa == [])
    return dict(view=view, giant=giant, alpha2=a2, shared_exon_index=S + 1,
                giant_exon127_len=up_e - up_s, first_intron_bp=intron_up, pairs=chosen)


def main():
    out, lines = {}, []
    lines.append("species\tgiant_transcript\tpair\tforward (5'->3')\tF_Tm\tF_GC%\tF_site\t"
                 "reverse (5'->3')\tR_Tm\tR_GC%\tR_site\tproduct_bp\tqPCR?\tin_silico")
    for view, (giant, a2) in SPECIES.items():
        r = design(view, giant, a2)
        out[view] = r
        print(f"\n### {LABEL[view]}  giant={giant}  (giant-specific exon 'e127' = {r['giant_exon127_len']} nt; "
              f"first intron {r['first_intron_bp']} bp)")
        for k, p in enumerate(r["pairs"], 1):
            f, rv = p["forward"], p["reverse"]
            fsite = "e127/shared junction" if f["junction"] else "giant exon 127 (giant-unique vs a2)"
            vflag = "OK" if p["validated"] else "CHECK"
            lines.append(f"{LABEL[view]}\t{giant}\t{k}\t{f['seq']}\t{f['tm']}\t{f['gc']}\t{fsite}\t"
                         f"{rv['seq']}\t{rv['tm']}\t{rv['gc']}\t{rv['loc']}\t{p['product']}\t"
                         f"{'yes' if p['qpcr_suitable'] else 'no'}\t{vflag}")
            print(f"  pair {k} [{vflag}]: F {f['seq']} (Tm {f['tm']}, {fsite})")
            print(f"            R {rv['seq']} (Tm {rv['tm']}, {rv['loc']})  -> {p['product']} bp"
                  f"{'  [qPCR]' if p['qpcr_suitable'] else ''}  giant={p['giant_products']} a2={p['alpha2_products']}")
    save_json(RESULTS / "mammalian_giant_primers.json", out)
    (RESULTS / "mammalian_giant_primers.tsv").write_text("\n".join(lines) + "\n")
    print(f"\nwrote {RESULTS/'mammalian_giant_primers.tsv'}")


if __name__ == "__main__":
    main()
