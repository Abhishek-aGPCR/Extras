#!/usr/bin/env python3
"""Step 11 - cross-species and developmental test of the nesprin-1alpha2 (muscle isoform) junction.

Why
  The adult zebrafish scan (step 10) is one tissue, one stage. Two checks make its result
  interpretable (both requested):
    1. Positive control - does the identical method detect the alpha2 junction where the isoform is
       known to be made? Human heart (every cardiomyocyte expresses nesprin-1alpha2; Holt 2019) and
       mouse C2C12 myotubes (alpha2 is switched on as myoblasts fuse; Stroud 2017).
    2. Development - alpha2 is developmentally dynamic (induced during myogenesis, then NMJ-restricted
       in adult skeletal muscle), so a time-restricted exon is easiest to see across differentiation
       / developmental stages, not in adult bulk muscle. Mouse myoblast vs myotube, and zebrafish
       8-12 hpf (early somites) vs 96 hpf / 4 dpf larvae (differentiated muscle).

What it does (identical method to step 10, made strand- and species-generic)
  For each panel: build capture k-mers over the alpha2/giant region of that species' locus, stream the
  FASTQ from ENA, keep reads sharing a k-mer, align them to the locus with minimap2 (preset splice:sr),
  and count - per sample and pooled per group - the alpha2 junction (exon 1alpha2 -> first shared exon)
  against the giant junction (previous giant exon -> the same shared exon), as independent fragments
  with duplicates removed, plus an alignment-free exact 24-mer count and exon 1alpha2 coverage. The
  giant junction is the built-in positive control (if it is detected and alpha2 is not, alpha2 is
  genuinely absent in that sample).

  Human SYNE1 and mouse Syne1 are on the minus strand; junctions are defined from the genomic intron
  between transcript-adjacent exons, so the counting is identical regardless of strand.

Panels (wild-type / control libraries only):
  human_heart   GSM9484163-165  healthy adult left ventricle (HC_LV 1-3)          [positive control]
  mouse_c2c12   PRJEB108758     C2C12 myoblast vs myotube day 3 (polysome RNA-seq)[dev switch + ctrl]
  zebrafish_dev PRJNA1483411 (8/10/12 hpf) + PRJEB98861 (96 hpf, 4 dpf WT larvae) [development]
  zebrafish_adult  reuses the step-10 captures (fast/slow muscle, ventricle) for the joint figure.

Usage
  python3 11_rnaseq_crossspecies.py patterns  [--panels ...]
  python3 11_rnaseq_crossspecies.py scan      [--panels ...] [--max-bytes N] [--workers 6]
  python3 11_rnaseq_crossspecies.py analyse   [--panels ...]
Requires curl, gzip, ripgrep (or grep) and mappy; shares helpers with 10_rnaseq_junctions.py.
"""
from __future__ import annotations

import argparse
import gzip
import importlib.util
import math
import re
import sys
import time
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from common import DATA, RESULTS, load_json, read_fasta, revcomp, save_json, ucsc_track  # noqa: E402

# reuse the tested download + alignment primitives from step 10
_spec = importlib.util.spec_from_file_location("rnaseq10", HERE / "10_rnaseq_junctions.py")
rj = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(rj)

K = rj.K_CAPTURE
RNADIR = DATA / "rnaseq"
OUT = RESULTS / "rnaseq"

# --- panels -------------------------------------------------------------------------------------
# samples: (run, label, group). Groups are pooled; the group order is the biological axis.
PANELS = {
    "human_heart": dict(
        view="human_SYNE1", giant="NM_182961.4", alpha2="NM_001347702.2",
        source="GEO GSM9484163-165 (healthy adult left ventricle)", kind="positive control",
        group_order=["left ventricle (healthy adult)"],
        samples=[("SRR37034261", "HC_LV_1", "left ventricle (healthy adult)"),
                 ("SRR37034259", "HC_LV_2", "left ventricle (healthy adult)"),
                 ("SRR37034257", "HC_LV_3", "left ventricle (healthy adult)")]),
    "human_brain": dict(
        view="human_SYNE1", giant="NM_182961.4", alpha2="NM_001347702.2",
        source="PRJNA1517665 (adult frontal cortex, QSBB; whole-tissue homogenate)",
        kind="non-muscle tissue (brain)",
        group_order=["frontal cortex"],
        patterns=RNADIR / "human_heart" / "capture_kmers.txt",   # same SYNE1 locus as human_heart
        samples=[("SRR40332934", "FCx_control", "frontal cortex"),
                 ("SRR40332935", "FCx_C9FTD_1", "frontal cortex"),
                 ("SRR40332938", "FCx_C9FTD_2", "frontal cortex"),
                 ("SRR40332939", "FCx_C9FTD_3", "frontal cortex")]),
    "human_osteoclast": dict(
        view="human_SYNE1", giant="NM_182961.4", alpha2="NM_001347702.2",
        source="PRJNA1504327 (monocyte-derived human osteoclasts)", kind="osteoclast",
        group_order=["osteoclast (monocyte-derived)"],
        patterns=RNADIR / "human_heart" / "capture_kmers.txt",   # same SYNE1 locus as human_heart
        samples=[("SRR39902634", "OC_mono_1", "osteoclast (monocyte-derived)"),
                 ("SRR39902635", "OC_mono_2", "osteoclast (monocyte-derived)"),
                 ("SRR39902636", "OC_mono_3", "osteoclast (monocyte-derived)")]),
    "mouse_c2c12": dict(
        view="mouse_Syne1", giant="NM_001347732.2", alpha2="NM_022027.3",
        source="PRJEB108758 (C2C12 myoblast vs myotube day 3, polysome RNA-seq)",
        kind="developmental switch + positive control",
        group_order=["myoblast", "myotube (day 3)"],
        samples=[("ERR16751628", "myoblast_rep1", "myoblast"),
                 ("ERR16749445", "myoblast_rep2", "myoblast"),
                 ("ERR16755182", "myotube_d3_rep1", "myotube (day 3)"),
                 ("ERR16755183", "myotube_d3_rep2", "myotube (day 3)")]),
    "zebrafish_dev": dict(
        view="zebrafish_syne1a", giant="ENSDART00000139350.4", alpha2="syne1a-alpha2-pred-plusDV23-GRCz11",
        source="PRJEB98861 (96 hpf / 4 dpf whole wild-type larvae; differentiated skeletal muscle)",
        kind="larval development",
        group_order=["96 hpf / 4 dpf larva"],
        samples=[("ERR15682390", "4dpf_wt_1", "96 hpf / 4 dpf larva"),
                 ("ERR15682391", "4dpf_wt_2", "96 hpf / 4 dpf larva"),
                 ("ERR15682392", "4dpf_wt_3", "96 hpf / 4 dpf larva"),
                 ("ERR15682393", "4dpf_wt_4", "96 hpf / 4 dpf larva")]),
    "zebrafish_larval_heart": dict(
        view="zebrafish_syne1a", giant="ENSDART00000139350.4", alpha2="syne1a-alpha2-pred-plusDV23-GRCz11",
        source="PRJNA395216 (dissected 48 hpf larval heart, wild-type/control siblings)",
        kind="larval heart",
        group_order=["48 hpf heart (control)"],
        patterns=RNADIR / "zebrafish_dev" / "capture_kmers.txt",   # identical syne1a locus as zebrafish_dev
        samples=[("SRR5855250", "48hpf_heart_ctrl_1", "48 hpf heart (control)"),
                 ("SRR5855251", "48hpf_heart_ctrl_2", "48 hpf heart (control)"),
                 ("SRR5855252", "48hpf_heart_ctrl_3", "48 hpf heart (control)"),
                 ("SRR5855253", "48hpf_heart_ctrl_4", "48 hpf heart (control)")]),
    "zebrafish_osteoclast": dict(
        view="zebrafish_syne1a", giant="ENSDART00000139350.4", alpha2="syne1a-alpha2-pred-plusDV23-GRCz11",
        source="PRJNA1331843 (Ctsk+ wild-type sorted osteoclasts)", kind="osteoclast",
        group_order=["osteoclast (Ctsk+ WT)"],
        samples=[("SRR35522660", "osteoclast_wt_1", "osteoclast (Ctsk+ WT)"),
                 ("SRR35522659", "osteoclast_wt_2", "osteoclast (Ctsk+ WT)"),
                 ("SRR35522658", "osteoclast_wt_3", "osteoclast (Ctsk+ WT)")]),
    # reuse step-10 captures (no re-download) so the joint figure spans all species/stages
    "zebrafish_adult": dict(
        view="zebrafish_syne1a", giant="ENSDART00000139350.4", alpha2="syne1a-alpha2-pred-plusDV23-GRCz11",
        source="PRJNA1257549 (adult fast/slow muscle, ventricle)", kind="adult tissue",
        group_order=["fast muscle", "slow muscle", "ventricle"],
        hits_dir=RNADIR / "hits", patterns=RNADIR / "capture_kmers.txt",
        samples=[(s[0], s[1], s[2]) for s in rj.SAMPLES]),
}


def panel_dir(name):
    return RNADIR / name


def hits_dir(name):
    return PANELS[name].get("hits_dir") or (panel_dir(name) / "hits")


def patterns_path(name):
    return PANELS[name].get("patterns") or (panel_dir(name) / "capture_kmers.txt")


# --- locus, model, key junctions (strand-generic) -----------------------------------------------
def intron_between(ex1, ex2):
    """Genomic intron between two transcript-adjacent exons, as (d, a) with d < a:
    d = last genomic base of the genomically-lower exon, a = first genomic base of the higher exon.
    This matches minimap2 N-gap coordinates regardless of gene strand."""
    if ex1[1] < ex2[0]:
        return (ex1[1], ex2[0])
    return (ex2[1], ex1[0])


def load_locus(name):
    P = PANELS[name]
    built = load_json(DATA / "models_built.json")[P["view"]]
    region = built["region"]
    seq = next(iter(read_fasta((DATA / "genomic" / f"{P['view']}.fa").read_text()).values())).upper()
    assert len(seq) == region["end"] - region["start"] + 1
    tx = {t["id"]: t for t in built["transcripts"]}
    giant = [tuple(e) for e in tx[P["giant"]]["exons"]]
    a2ex = [tuple(e) for e in tx[P["alpha2"]]["exons"]]
    strand = tx[P["giant"]]["strand"]
    shared = a2ex[1]                                   # first exon shared between alpha2 and giant
    assert shared in giant, f"{name}: alpha2 second exon is not a giant exon"
    i = giant.index(shared)
    key = {
        "alpha2": dict(intron=intron_between(a2ex[0], shared), label="1α2>shared"),
        "giant": dict(intron=intron_between(giant[i - 1], shared), label="giant>shared"),
    }
    other = {k: [tuple(e) for e in v["exons"]] for k, v in tx.items()}
    return dict(name=name, view=P["view"], region=region, seq=seq, start=region["start"], strand=strand,
                giant=giant, alpha2_ex=a2ex, exon1a2=a2ex[0], shared=shared, shared_idx=i,
                giant_prev=giant[i - 1], key=key, other=other,
                ucsc_db=built["view"]["ucsc_db"], chrom=built["view"]["chrom"])


def norm_junc(seq, start, d, a):
    """Normalise a genomic junction (d<a) to the preferred motif position; return (d,a,motif) genomic."""
    d2, a2, motif = rj.shift_junction(seq, d - start, a - start)
    return d2 + start, a2 + start, motif


def tx_motif(seq, start, intron, strand):
    """The true splice motif in transcription orientation (e.g. GT-AG, GC-AG) from the annotated intron
    boundaries, so minus-strand genes read correctly rather than as their + strand complement."""
    d, a = intron                                   # d<a genomic; intron bases d+1..a-1
    donor = seq[d + 1 - start:d + 3 - start]
    acc = seq[a - 2 - start:a - start]
    return (donor + "-" + acc) if strand == 1 else (revcomp(acc) + "-" + revcomp(donor))


def features(L):
    start, seq = L["start"], L["seq"]
    key = {}
    for k, v in L["key"].items():
        d, a = v["intron"]
        dn, an, m = norm_junc(seq, start, d, a)
        key[k] = dict(donor=dn, acceptor=an, motif=m, motif_tx=tx_motif(seq, start, (d, a), L["strand"]),
                      intron=(d, a))
    # annotated genomic introns (any transcript) for the anchor rule
    annotated = set()
    for ex in L["other"].values():
        for i in range(len(ex) - 1):
            d, a = intron_between(ex[i], ex[i + 1])
            annotated.add((d, a))
    for v in key.values():
        annotated.add((v["donor"], v["acceptor"]))
    # capture intervals: giant exons near the shared exon, exon 1alpha2, and a window around it
    i = L["shared_idx"]
    neigh = L["giant"][max(0, i - 3): i + 4]
    a2 = L["exon1a2"]
    win = (min(a2[0], L["shared"][0]) - 1500, max(a2[1], L["shared"][1]) + 1500)
    intervals = merge_intervals(neigh + [a2, L["alpha2_ex"][1]] + [win])
    exonic = set()
    for s, e in neigh + list(L["alpha2_ex"][:3]):
        exonic.update(range(s, e + 1))
    return dict(key=key, annotated=annotated, intervals=intervals, exonic=exonic, neigh=neigh)


def merge_intervals(iv):
    iv = sorted((min(a, b), max(a, b)) for a, b in iv)
    out = [list(iv[0])]
    for a, b in iv[1:]:
        if a <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return [tuple(x) for x in out]


# --- patterns -----------------------------------------------------------------------------------
def make_patterns(name):
    L = load_locus(name)
    F = features(L)
    start, seq = L["start"], L["seq"]
    lo = min(a for a, _ in F["intervals"]) - 1
    hi = max(b for _, b in F["intervals"])
    rm = ucsc_track(L["ucsc_db"], "rmsk", "chr" + str(L["chrom"]), lo, hi)
    masked = [(r["genoStart"] + 1, r["genoEnd"]) for r in rm]
    exonic = F["exonic"]

    def is_masked(a, b):
        return any(not (b < m0 or a > m1) for m0, m1 in masked) and not all(x in exonic for x in range(a, b + 1))

    kmers, drep, dlc = set(), 0, 0
    for a, b in F["intervals"]:
        for g in range(a, b - K + 2):
            k = seq[g - start: g - start + K]
            if "N" in k:
                continue
            if is_masked(g, g + K - 1):
                drep += 1
                continue
            if rj._entropy(k) < 1.5 or re.search(r"A{10}|C{10}|G{10}|T{10}", k):
                dlc += 1
                continue
            kmers.add(k)
            kmers.add(revcomp(k))
    panel_dir(name).mkdir(parents=True, exist_ok=True)
    patterns_path(name).write_text("\n".join(sorted(kmers)) + "\n")
    save_json(panel_dir(name) / "capture_kmers.json",
              dict(k=K, n_patterns=len(kmers), intervals=F["intervals"], repeats_masked=masked,
                   dropped_repeat=drep, dropped_low_complexity=dlc))
    print(f"[{name}] capture k-mers: {len(kmers)} over {len(F['intervals'])} intervals "
          f"(dropped {drep} repeat, {dlc} low-complexity)")


# --- scan (reuses step-10 download + filter) ----------------------------------------------------
def scan(name, max_bytes, workers, chunk, mates=(1, 2)):
    if not patterns_path(name).exists():
        make_patterns(name)
    hd = hits_dir(name)
    hd.mkdir(parents=True, exist_ok=True)
    rj.PATTERNS = patterns_path(name)       # process_file / grep read these globals
    rj.HITS = hd
    todo = []
    for run, label, group in PANELS[name]["samples"]:
        for mate in mates:
            info = hd / f"{run}_{mate}.json"
            if info.exists():
                old = load_json(info)
                if old.get("complete") or (max_bytes and old.get("bytes_read", 0) >= max_bytes):
                    continue
            todo.append((run, mate))
    with ThreadPoolExecutor(workers) as dl:
        for run, mate in todo:
            url = rj.ena_url(run, mate)
            try:
                size = rj.remote_size(url)
            except Exception as e:
                print(f"[{name}] {run}_{mate}: no file at ENA ({e}); skipping (single-end?)")
                continue
            nbytes = min(size, max_bytes) if max_bytes else size
            d = rj.PARTS / f"{run}_{mate}"
            d.mkdir(parents=True, exist_ok=True)
            ranges = [(a, min(a + chunk, nbytes) - 1) for a in range(0, nbytes, chunk)]
            parts = [d / f"part_{i:05d}" for i in range(len(ranges))]
            t0 = time.time()
            for f in [dl.submit(rj.fetch_range, url, a, b, p) for (a, b), p in zip(ranges, parts)]:
                f.result()
            print(f"[{name}] {run}_{mate}: {nbytes / 1e6:,.0f} MB in {(time.time() - t0) / 60:.1f} min", flush=True)
            try:
                rj.process_file(run, mate, parts, nbytes == size, size, nbytes)
            except Exception as e:
                print(f"[{name}] {run}_{mate}: processing failed ({str(e)[:120]}); skipping this file", flush=True)
                for p in parts:
                    p.unlink(missing_ok=True)


def read_hits(name, run):
    hd = hits_dir(name)
    for mate in (1, 2):
        p = hd / f"{run}_{mate}.tsv.gz"
        if not p.exists():
            continue
        with gzip.open(p, "rt") as fh:
            for line in fh:
                nm, seq = line.rstrip("\n").split("\t")
                frag = nm.lstrip("@").split()[0]
                if frag.endswith("/1") or frag.endswith("/2"):
                    frag = frag[:-2]
                yield mate, frag, seq.upper()


# --- analyse ------------------------------------------------------------------------------------
def analyse_panel(name, al_cache):
    import mappy
    L = load_locus(name)
    F = features(L)
    start, seq = L["start"], L["seq"]
    al = al_cache.get(L["view"])
    if al is None:
        al = mappy.Aligner(seq=seq, preset="splice:sr")
        al_cache[L["view"]] = al
    key = F["key"]
    a2key = (key["alpha2"]["donor"], key["alpha2"]["acceptor"])
    gkey = (key["giant"]["donor"], key["giant"]["acceptor"])
    # exact 24-mers for the two key junctions
    j24 = {}
    for k, v in key.items():
        d, a = v["donor"] - start, v["acceptor"] - start
        s = seq[d - 11:d + 1] + seq[a:a + 12]
        j24[k] = (s, revcomp(s))
    a2s, a2e = L["exon1a2"][0] - start, L["exon1a2"][1] - start
    sh_s, sh_e = L["shared"][0] - start, L["shared"][1] - start
    resc_cands = {"alpha2": a2key[0] if a2key[0] == gkey[0] else a2key[1],
                  "giant": gkey[0] if a2key[0] == gkey[0] else gkey[1]}
    shared_boundary = (set(a2key) & set(gkey))
    rescue_ok = len(shared_boundary) == 1
    present = [s for s in PANELS[name]["samples"] if any((hits_dir(name) / f"{s[0]}_{m}.json").exists() for m in (1, 2))]
    rows = []
    for run, label, group in present:
        infos = [load_json(hits_dir(name) / f"{run}_{m}.json") for m in (1, 2)
                 if (hits_dir(name) / f"{run}_{m}.json").exists()]
        frag_j, frag_ju, frag24 = defaultdict(set), defaultdict(set), defaultdict(set)
        seen_aln, seen_seq = set(), set()
        n_cap = n_map = 0
        a2cov = Counter()
        into_a2 = set()
        for mate, frag, rd in read_hits(name, run):
            n_cap += 1
            sseq = (mate, min(rd, revcomp(rd)))
            dupseq = sseq in seen_seq
            seen_seq.add(sseq)
            for k, (s1, s2) in j24.items():
                if s1 in rd or s2 in rd:
                    if not dupseq:
                        frag24[k].add(frag)
            r = rj.align_read(al, seq, rd)
            if not r or not r["ok"]:
                continue
            n_map += 1
            sig = (mate, r["strand"], r["blocks"][0][0], r["blocks"][-1][1], tuple(j[:2] for j in r["juncs"]))
            dup = sig in seen_aln
            seen_aln.add(sig)
            for b0, b1 in r["blocks"]:
                o = min(b1, a2e + 1) - max(b0, a2s)
                if o > 0:
                    a2cov["exon1a2"] += o
                o = min(b1, sh_e + 1) - max(b0, sh_s)
                if o > 0:
                    a2cov["shared"] += o
            for d, a, motif, left, right, near, qj in r["juncs"]:
                g = (d + start, a + start)
                need = rj.MIN_ANCHOR if g in F["annotated"] or g in (a2key, gkey) else rj.MIN_ANCHOR_NOVEL
                if min(left, right) < need or near > 1:
                    continue
                frag_j[g].add(frag)
                if not dup:
                    frag_ju[g].add(frag)
                if (a2s <= d <= a2e or a2s <= a <= a2e) and g != a2key:
                    into_a2.add(frag)
            if rescue_ok:
                rk = rj.rescue_clipped(r, seq, resc_cands, list(shared_boundary)[0] - start)
                if rk:
                    g = a2key if rk == "alpha2" else gkey
                    frag_j[g].add(frag)
                    if not dup:
                        frag_ju[g].add(frag)
        n_reads = sum(i["reads_scanned"] for i in infos)
        rows.append(dict(
            panel=name, run=run, label=label, group=group, mates=len(infos),
            complete=len(infos) == max(len(infos), 2) and all(i["complete"] for i in infos),
            reads_scanned=n_reads, reads_captured=n_cap, reads_mapped=n_map,
            alpha2_junction=len(frag_ju.get(a2key, ())), alpha2_junction_24mer=len(frag24["alpha2"]),
            giant_junction=len(frag_ju.get(gkey, ())), giant_junction_24mer=len(frag24["giant"]),
            junctions_into_exon1a2=len(into_a2),
            exon1a2_mean_cov=round(a2cov["exon1a2"] / max(1, a2e - a2s + 1), 2),
            shared_exon_mean_cov=round(a2cov["shared"] / max(1, sh_e - sh_s + 1), 2)))
        print(f"[{name}] {run} {label} ({group}): captured {n_cap:,} mapped {n_map:,}; "
              f"alpha2 {rows[-1]['alpha2_junction']} (24-mer {rows[-1]['alpha2_junction_24mer']}) "
              f"giant {rows[-1]['giant_junction']}; exon1a2 cov {rows[-1]['exon1a2_mean_cov']} "
              f"shared cov {rows[-1]['shared_exon_mean_cov']}", flush=True)
    groups = []
    for grp in PANELS[name]["group_order"]:
        gr = [r for r in rows if r["group"] == grp]
        if not gr:
            continue
        a = sum(r["alpha2_junction"] for r in gr)
        g = sum(r["giant_junction"] for r in gr)
        nreads = sum(r["reads_scanned"] for r in gr)
        groups.append(dict(
            group=grp, n_samples=len(gr), reads_scanned=nreads,
            alpha2_junction=a, alpha2_junction_24mer=sum(r["alpha2_junction_24mer"] for r in gr),
            giant_junction=g, giant_junction_24mer=sum(r["giant_junction_24mer"] for r in gr),
            alpha2_fraction=round(a / (a + g), 4) if a + g else None,
            alpha2_fraction_ci95=wilson(a, a + g),
            alpha2_per_10M=round(a / nreads * 1e7, 3) if nreads else None,
            giant_per_10M=round(g / nreads * 1e7, 3) if nreads else None,
            exon1a2_mean_cov=round(sum(r["exon1a2_mean_cov"] for r in gr) / len(gr), 2),
            shared_exon_mean_cov=round(sum(r["shared_exon_mean_cov"] for r in gr) / len(gr), 2),
            junctions_into_exon1a2=sum(r["junctions_into_exon1a2"] for r in gr)))
    P = PANELS[name]
    meta = dict(source=P["source"], kind=P["kind"], view=P["view"], group_order=P["group_order"])
    return dict(panel=name, meta=meta, view=L["view"], strand=L["strand"],
                coords=dict(assembly=load_json(DATA / "models_built.json")[L["view"]]["view"]["assembly"],
                            chrom=L["chrom"], exon1a2=list(L["exon1a2"]), shared=list(L["shared"]),
                            alpha2_intron=list(key["alpha2"]["intron"]), alpha2_motif=key["alpha2"]["motif_tx"],
                            giant_motif=key["giant"]["motif_tx"]),
                samples=rows, groups=groups)


def wilson(k, n, z=1.96):
    if not n:
        return [None, None]
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return [round(max(0.0, c - h), 4), round(min(1.0, c + h), 4)]


def analyse(names):
    OUT.mkdir(parents=True, exist_ok=True)
    al_cache = {}
    results = {}
    for name in names:
        if not any((hits_dir(name) / f"{s[0]}_1.json").exists() for s in PANELS[name]["samples"]):
            print(f"[{name}] no scanned reads; run scan first — skipping")
            continue
        results[name] = analyse_panel(name, al_cache)
    if not results:
        raise SystemExit("nothing analysed")
    save_json(OUT / "crossspecies_summary.json", results)
    # flat TSV
    with open(OUT / "crossspecies_groups.tsv", "w") as fh:
        cols = ["panel", "kind", "group", "n_samples", "reads_scanned", "alpha2_junction", "alpha2_junction_24mer",
                "giant_junction", "alpha2_fraction", "alpha2_per_10M", "giant_per_10M",
                "exon1a2_mean_cov", "shared_exon_mean_cov", "junctions_into_exon1a2"]
        fh.write("\t".join(cols) + "\n")
        for name, R in results.items():
            for g in R["groups"]:
                fh.write("\t".join(str({**g, "panel": name, "kind": R["meta"]["kind"]}.get(c, "")) for c in cols) + "\n")
    (OUT / "crossspecies.svg").write_text(figure(results))
    print("\n=== cross-species / developmental summary ===")
    for name, R in results.items():
        print(f"\n{name} [{R['meta']['kind']}] — {R['meta']['source']}")
        print(f"  alpha2 intron {R['coords']['exon1a2']} motif {R['coords']['alpha2_motif']}; {R['strand']:+d} strand")
        for g in R["groups"]:
            print(f"    {g['group']:32s} alpha2={g['alpha2_junction']:4d} (24mer {g['alpha2_junction_24mer']:4d}) "
                  f"giant={g['giant_junction']:5d} share={g['alpha2_fraction']} "
                  f"alpha2/10M={g['alpha2_per_10M']} exon1a2cov={g['exon1a2_mean_cov']} sharedcov={g['shared_exon_mean_cov']}")


# --- figure -------------------------------------------------------------------------------------
def figure(results):
    order = [n for n in ["human_heart", "human_brain", "human_osteoclast", "mouse_c2c12",
                         "zebrafish_larval_heart", "zebrafish_dev", "zebrafish_adult", "zebrafish_osteoclast"] if n in results]
    rowsflat = []
    for name in order:
        R = results[name]
        for g in R["groups"]:
            rowsflat.append((name, R, g))
    W, left, right, top, rowh = 1280, 290, 410, 70, 30
    H = top + len(rowsflat) * rowh + 80
    maxlog = max([math.log10(1 + g["giant_per_10M"] + g["alpha2_per_10M"]) for _, _, g in rowsflat] + [1])

    def bar(x, w, col):
        return f'<rect x="{x:.1f}" y="0" width="{max(w,0.4):.1f}" height="12" fill="{col}"/>'
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" '
         'font-family="Helvetica,Arial,sans-serif" font-size="12">',
         f'<rect width="{W}" height="{H}" fill="#ffffff"/>',
         '<text x="24" y="30" font-size="16" font-weight="bold">Nesprin-1α2 junction across species, tissues and development</text>',
         '<text x="24" y="50" fill="#555" font-size="12">Junction-spanning fragments per 10M reads (log scale): '
         'red = α2 (exon 1α2→shared exon), grey = giant. Counts and α2 share at right.</text>']
    y = top
    scale = (W - left - right) / maxlog
    for name, R, g in rowsflat:
        lab = f"{R['meta']['view'].split('_')[0]} · {g['group']}"
        o.append(f'<text x="24" y="{y + 10}" font-size="11.5">{lab[:46]}</text>')
        gx = math.log10(1 + g["giant_per_10M"]) * scale
        ax = math.log10(1 + g["alpha2_per_10M"]) * scale
        o.append(f'<g transform="translate({left},{y})">' + bar(0, gx, "#8a95a3") + '</g>')
        o.append(f'<g transform="translate({left},{y})">' + bar(0, ax, "#b0124c") + '</g>')
        frac = '–' if g["alpha2_fraction"] is None else f'{g["alpha2_fraction"]*100:.2f}%'
        detect = "α2 detected" if g["alpha2_junction"] >= 3 else ("trace" if g["alpha2_junction"] else "not detected")
        col = "#1d7a4c" if g["alpha2_junction"] >= 3 else ("#a15c00" if g["alpha2_junction"] else "#8a95a3")
        o.append(f'<text x="{left + (W - left - right) + 10}" y="{y + 10}" font-size="11.5" fill="{col}" font-weight="600">'
                 f'{detect} · α2 {g["alpha2_junction"]} vs giant {g["giant_junction"]} · {frac}</text>')
        y += rowh
    # panel separators + labels
    yy = top
    for name in order:
        n = len(results[name]["groups"])
        o.append(f'<line x1="24" x2="{W-24}" y1="{yy-6}" y2="{yy-6}" stroke="#e5e8ec"/>')
        o.append(f'<text x="24" y="{yy-10}" font-size="10.5" fill="#888">{results[name]["meta"]["kind"]}</text>')
        yy += n * rowh
    o.append(f'<text x="24" y="{H-20}" fill="#888" font-size="10">Method: capture + minimap2 splice:sr, '
             'fragments with duplicates removed; giant junction is the internal positive control. '
             'Predictions from public RNA-seq; see results/rnaseq/crossspecies_groups.tsv.</text>')
    o.append("</svg>")
    return "\n".join(o)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for c in ("patterns", "scan", "analyse"):
        s = sub.add_parser(c)
        s.add_argument("--panels", nargs="*", default=["human_heart", "mouse_c2c12", "zebrafish_dev"])
        if c == "scan":
            s.add_argument("--max-bytes", type=int, default=None)
            s.add_argument("--workers", type=int, default=6)
            s.add_argument("--chunk", type=int, default=32 * 1024 * 1024)
    a = ap.parse_args()
    if a.cmd == "patterns":
        for n in a.panels:
            make_patterns(n)
    elif a.cmd == "scan":
        for n in a.panels:
            scan(n, a.max_bytes, a.workers, a.chunk)
    elif a.cmd == "analyse":
        analyse(a.panels)


if __name__ == "__main__":
    main()
