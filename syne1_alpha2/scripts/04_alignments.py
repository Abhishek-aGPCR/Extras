#!/usr/bin/env python3
"""Step 4 - precomputed alignments for the viewer (FASTA copies in results/alignments/).

  alpha2_protein        nesprin-1 alpha2 proteins: human, mouse, rat, zebrafish (predicted) +
                        representative bird, amphibian and fish alpha2 RefSeq models (FAMSA)
  alpha2_cds_codon      codon alignment of alpha2 CDS (human/mouse/rat/zebrafish), from the protein MSA
  exon1alpha2_nt        exon 1alpha2 + 12 nt of the downstream intron (GC donor), center-star alignment
  alpha2_nterm_vertebrates  unique N-termini of all alpha2-type RefSeq proteins (survey) + zebrafish
  giant_protein         giant isoforms of human, mouse, rat, zebrafish syne1a and syne1b (FAMSA)
  alpha2_exon_table     exon-by-exon correspondence of the alpha2 isoforms with per-exon identity
"""
from __future__ import annotations

import sys
from pathlib import Path

import parasail
from pyfamsa import Aligner, Sequence

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DATA, RESULTS, load_json, ncbi_efetch, read_fasta, revcomp, save_json, write_fasta  # noqa: E402

ALN_DIR = RESULTS / "alignments"
ALN_DIR.mkdir(exist_ok=True)
BUILT = load_json(DATA / "models_built.json")
ANA = load_json(RESULTS / "alpha2_exon_analysis.json")


def tx(key, tid):
    for t in BUILT[key]["transcripts"]:
        if t["id"].split(".")[0] == tid.split(".")[0] or t["id"] == tid:
            return t
    raise KeyError(tid)


def genome(key):
    return "".join(open(DATA / "genomic" / f"{key}.fa").read().split("\n")[1:])


def sub(key, a, b, strand):
    st = BUILT[key]["region"]["start"]
    s = genome(key)[a - st:b - st + 1]
    return s if strand == 1 else revcomp(s)


def cds_of(key, t):
    a, b = t["cds"]
    return "".join(sub(key, max(s, a), min(e, b), t["strand"]) for s, e in t["exons"] if max(s, a) <= min(e, b))


def famsa(records):
    al = Aligner(guide_tree="upgma")
    msa = al.align([Sequence(n.encode(), s.encode()) for n, s in records])
    order = {n: i for i, (n, _) in enumerate(records)}
    rows = sorted(((s.id.decode(), s.sequence.decode()) for s in msa), key=lambda r: order[r[0]])
    return rows


def center_star(center, others, match=5, mismatch=-4, go=10, ge=1):
    """Nucleotide MSA by semi-global alignment of every sequence to a center sequence."""
    mat = parasail.matrix_create("ACGTN", match, mismatch)
    cname, cseq = center
    pair = []
    for name, s in others:
        r = parasail.sg_trace_scan_16(s, cseq, go, ge, mat)
        pair.append((name, r.traceback.query, r.traceback.ref))
    # insertion lengths before each center position (0..len)
    ins = [0] * (len(cseq) + 1)
    per = []
    for name, q, ref in pair:
        pos, cur, d = 0, 0, {}
        for a, b in zip(q, ref):
            if b == "-":
                cur += 1
            else:
                d[pos] = cur
                ins[pos] = max(ins[pos], cur)
                cur = 0
                pos += 1
        d[pos] = cur
        ins[pos] = max(ins[pos], cur)
        per.append((name, q, ref))
    # build rows
    def expand(q, ref):
        out, pos, buf = [], 0, []
        for a, b in zip(q, ref):
            if b == "-":
                buf.append(a)
            else:
                out.append("".join(buf).ljust(ins[pos], "-"))
                out.append(a)
                buf = []
                pos += 1
        out.append("".join(buf).ljust(ins[pos], "-"))
        return "".join(out)
    rows = [(cname, expand(cseq, cseq))]
    for name, q, ref in per:
        rows.append((name, expand(q, ref)))
    return rows


def save(name, rows, meta=None):
    write_fasta(ALN_DIR / f"{name}.fa", rows)
    return dict(name=name, rows=[[n, s] for n, s in rows], **(meta or {}))


def main() -> None:
    out = {}
    # ---------------- alpha2 proteins ----------------
    core = [("Human_NP_001334631.1", tx("human_SYNE1", "NM_001347702.2")["protein"]),
            ("Mouse_XP_036011830.1", tx("mouse_Syne1", "XM_036155937.1")["protein"]),
            ("Rat_XP_063126704.1", tx("rat_Syne1", "XM_063270634.2")["protein"]),
            ("Zebrafish_syne1a_alpha2_pred", tx("zebrafish_syne1a", "syne1a-alpha2-pred-plusDV23-GRCz11")["protein"])]
    extra_ids = {"Chicken_XP_040522257.2": "XP_040522257.2", "Xenopus_XP_031757385.1": "XP_031757385.1",
                 "Herring_XP_042565831.1": "XP_042565831.1", "Pike_syne1a_XP_028970453.2": "XP_028970453.2",
                 "Betta_syne1a_XP_055362111.1": "XP_055362111.1", "Sterlet_syne1a_XP_033873712.3": "XP_033873712.3"}
    fa = read_fasta(ncbi_efetch("protein", list(extra_ids.values())))
    recs = core + [(k, fa[v]) for k, v in extra_ids.items()]
    rows = famsa(recs)
    out["alpha2_protein"] = save("alpha2_protein", rows, dict(
        type="protein", title="Nesprin-1α2 proteins (human, mouse, rat, zebrafish predicted + other vertebrates)"))

    # ---------------- alpha2 CDS codon alignment ----------------
    cds = {"Human_NP_001334631.1": cds_of("human_SYNE1", tx("human_SYNE1", "NM_001347702.2")),
           "Mouse_XP_036011830.1": cds_of("mouse_Syne1", tx("mouse_Syne1", "XM_036155937.1")),
           "Rat_XP_063126704.1": cds_of("rat_Syne1", tx("rat_Syne1", "XM_063270634.2")),
           "Zebrafish_syne1a_alpha2_pred": cds_of("zebrafish_syne1a", tx("zebrafish_syne1a", "syne1a-alpha2-pred-plusDV23-GRCz11"))}
    codon_rows = []
    for name, aln in rows:
        if name not in cds:
            continue
        nt, i, s = cds[name], 0, []
        for c in aln:
            if c == "-":
                s.append("---")
            else:
                s.append(nt[3 * i:3 * i + 3])
                i += 1
        codon_rows.append((name, "".join(s)))
    out["alpha2_cds_codon"] = save("alpha2_cds_codon", codon_rows, dict(
        type="codon", title="Nesprin-1α2 coding sequences, codon alignment (from protein MSA)"))

    # ---------------- exon 1alpha2 nucleotide alignment ----------------
    mam = ANA["mammal_exon1alpha2"]
    zx = ANA["zebrafish_exon1alpha2"]

    def ex_with_donor(key, e, strand):
        if strand == -1:
            return sub(key, e[0] - 12, e[1], -1)
        return sub(key, e[0], e[1] + 12, 1)
    ntrecs = [("Human_exon1a2", ex_with_donor("human_SYNE1", mam["human"]["exon"], -1)),
              ("Mouse_exon1a2", ex_with_donor("mouse_Syne1", mam["mouse"]["exon"], -1)),
              ("Rat_exon1a2", ex_with_donor("rat_Syne1", mam["rat"]["exon"], -1)),
              ("Zebrafish_syne1a_exon1a2", ex_with_donor("zebrafish_syne1a", zx["predicted_exon"], 1))]
    nrows = center_star(ntrecs[0], ntrecs[1:])
    # annotation: ATG start and donor for each row
    ann = {}
    for (name, raw), (_, aln) in zip(ntrecs, nrows):
        sp = name.split("_")[0].lower()
        utr = mam[sp]["utr5"] if sp in mam else len(zx["exon_seq"]) - zx["coding_nt"]
        exon_len = len(raw) - 12
        cols, i = {}, 0
        for c_i, c in enumerate(aln):
            if c != "-":
                cols[i] = c_i
                i += 1
        ann[name] = dict(atg_col=cols.get(utr), donor_col=cols.get(exon_len), cds_end_col=cols.get(exon_len - 1))
    out["exon1alpha2_nt"] = save("exon1alpha2_nt", nrows, dict(
        type="dna", title="Exon 1α2 (5'UTR + coding part) and the first 12 nt of the downstream intron",
        annotations=ann))

    # ---------------- vertebrate unique N-termini ----------------
    surv = ANA["vertebrate_survey"]
    vrecs = [("Danio_rerio_syne1a_PREDICTED", zx["peptide"] + "R")]
    seen = set()
    for r in surv:
        nm = f"{r['organism'].replace(' ', '_')}_{r['accession']}"
        if nm in seen:
            continue
        seen.add(nm)
        vrecs.append((nm, r["unique_nterm"] + r["junction"][:1]))
    vrows = famsa(vrecs)
    groups = {f"{r['organism'].replace(' ', '_')}_{r['accession']}": r["group"] for r in surv}
    groups["Danio_rerio_syne1a_PREDICTED"] = "Ray-finned fish"
    order = ["Ray-finned fish", "Amphibians", "Reptiles", "Birds", "Mammals", "Other vertebrates"]
    vrows.sort(key=lambda r: (order.index(groups.get(r[0], "Other vertebrates")),
                              0 if "PREDICTED" in r[0] else 1, r[0]))
    out["alpha2_nterm_vertebrates"] = save("alpha2_nterm_vertebrates", vrows, dict(
        type="protein", title=f"Unique α2 N-termini of {len(vrows) - 1} vertebrate RefSeq proteins + zebrafish prediction "
                              "(last residue = codon split with the first shared exon)",
        groups={n: groups.get(n, "Other vertebrates") for n, _ in vrows}))

    # ---------------- giant proteins ----------------
    grecs = [("Human_NP_892006.3", tx("human_SYNE1", "NM_182961.4")["protein"]),
             ("Mouse_NP_001334661.1", tx("mouse_Syne1", "NM_001347732.2")["protein"]),
             ("Rat_NP_001406778.1", tx("rat_Syne1", "NM_001419849.1")["protein"]),
             ("Zebrafish_syne1a_ENSDARP_syne1a-201", tx("zebrafish_syne1a", "ENSDART00000139350")["protein"]),
             ("Zebrafish_syne1b_ENSDARP_syne1b-202", tx("zebrafish_syne1b", "ENSDART00000148520")["protein"])]
    out["giant_protein"] = save("giant_protein", famsa(grecs), dict(
        type="protein", title="Nesprin-1 giant isoforms (human, mouse, rat, zebrafish syne1a & syne1b)"))

    # ---------------- alpha2 exon correspondence ----------------
    tmap = {"Human_NP_001334631.1": ("human_SYNE1", "NM_001347702.2"),
            "Mouse_XP_036011830.1": ("mouse_Syne1", "XM_036155937.1"),
            "Rat_XP_063126704.1": ("rat_Syne1", "XM_063270634.2"),
            "Zebrafish_syne1a_alpha2_pred": ("zebrafish_syne1a", "syne1a-alpha2-pred-plusDV23-GRCz11")}
    prow = {n: s for n, s in rows}
    col_of = {}
    for n in tmap:
        c, i = [], 0
        for k, ch in enumerate(prow[n]):
            if ch != "-":
                c.append(k)
        col_of[n] = c
    table = []
    ht = tx(*tmap["Human_NP_001334631.1"])
    for xi, (s, e) in enumerate(ht["exons"]):
        aa = ht["exon_aa"][xi]
        row = dict(exon=xi + 1, human=dict(coords=[s, e], len=e - s + 1, aa=aa))
        if aa:
            hcols = col_of["Human_NP_001334631.1"]
            aa = [aa[0], min(aa[1], len(hcols)), aa[2]]
            c0, c1 = hcols[aa[0] - 1], hcols[aa[1] - 1]
            for n, (key, tid) in tmap.items():
                if n.startswith("Human"):
                    continue
                t = tx(key, tid)
                # the exon of t sharing the most alignment columns with this human exon
                best, best_ov = None, 0
                for yi, ya in enumerate(t["exon_aa"]):
                    if not ya:
                        continue
                    y1 = min(ya[1], len(col_of[n]))  # last exon's range includes the stop codon
                    if ya[0] > y1:
                        continue
                    d0, d1 = col_of[n][ya[0] - 1], col_of[n][y1 - 1]
                    ov = min(d1, c1) - max(d0, c0) + 1
                    if ov > best_ov:
                        best, best_ov = yi, ov
                ident = sum(1 for k in range(c0, c1 + 1) if prow[n][k] == prow["Human_NP_001334631.1"][k] != "-")
                a_, b_ = t["exons"][best]
                row[n.split("_")[0].lower()] = dict(exon=best + 1, coords=[a_, b_], len=b_ - a_ + 1,
                                                    identity=round(100 * ident / (aa[1] - aa[0] + 1), 1))
        table.append(row)
    out["alpha2_exon_table"] = table
    save_json(DATA / "alignments.json", out)
    for k, v in out.items():
        if isinstance(v, dict) and "rows" in v:
            print(f"{k}: {len(v['rows'])} sequences x {len(v['rows'][0][1])} columns")
    print("exon table rows:", len(table))
    for r in table:
        print(r["exon"], r["human"]["len"], {k: (v.get("exon"), v.get("len"), v.get("identity")) for k, v in r.items()
                                              if k not in ("exon", "human")})


if __name__ == "__main__":
    main()
