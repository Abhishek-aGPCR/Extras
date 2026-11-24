#!/usr/bin/env python3
"""Step 10 - test the predicted zebrafish exon 1alpha2 (syne1a) with public adult RNA-seq.

Question
  Is the exon predicted by homology in step 2 (syne1a, GRCz11 chr20:26,341,126-26,341,308, + strand)
  transcribed and spliced, through its GC-AG donor at 26,341,308, onto giant exon 128
  (26,342,757) in adult zebrafish muscle - i.e. does zebrafish make a nesprin-1alpha2-type mRNA?

Data
  BioProject PRJNA1257549: poly(A) RNA-seq, 2 x 150 bp, adult zebrafish fast-twitch muscle,
  slow-twitch muscle and cardiac ventricle, wild type (wt) and mb knockout (ko). Only the
  wild-type libraries are used. FASTQ files are read from ENA (ftp.sra.ebi.ac.uk, HTTPS).

Method (targeted; no genome-wide alignment)
  patterns  every 20-mer of syne1a exons 126-144 (giant model ENSDART00000139350.4) and of the
            whole of intron 127 (which hosts exon 1alpha2), minus RepeatMasker and low-complexity
            k-mers, both strands -> data/rnaseq/capture_kmers.txt
  scan      each FASTQ is downloaded as parallel byte ranges, decompressed, and every read that
            contains a capture k-mer is kept (grep -F) -> data/rnaseq/hits/<run>_<mate>.tsv.gz
            (+ .json with the number of reads scanned and whether the file was read completely).
  analyse   captured reads are aligned to the syne1a locus (GRCz11 chr20:26,093,467-26,361,438)
            with minimap2 (mappy, preset splice:sr). Junctions are taken from the CIGAR, shifted to
            the preferred motif among equivalent positions (GT-AG > GC-AG > AT-AC), and counted once
            per fragment (both mates of a fragment count once) with >= 8 nt aligned on each side.
            An alignment-free count (exact 24-mer spanning the junction, 12 + 12 nt) is reported
            for the key junctions as an independent check. Coverage, DV23 inclusion and the
            5' extent of exon 1alpha2 are also measured.
            -> results/rnaseq/*.tsv, rnaseq_summary.json, alpha2_junction_examples.txt,
               syne1a_alpha2_rnaseq.svg

Usage
  python3 10_rnaseq_junctions.py patterns
  python3 10_rnaseq_junctions.py scan [--runs SRR...] [--max-bytes N] [--workers 6]
  python3 10_rnaseq_junctions.py analyse
Requires: curl, gzip, mawk (or awk), ripgrep (or GNU grep), and the Python package mappy (pip install mappy).
"""
from __future__ import annotations

import argparse
import gzip
import math
import re
import shutil
import subprocess
import sys
import time
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DATA, EUTILS, RESULTS, http, load_json, read_fasta, revcomp, save_json, ucsc_track  # noqa: E402

VIEW = "zebrafish_syne1a"
GIANT = "ENSDART00000139350.4"
ALPHA2_PLUS = "syne1a-alpha2-pred-plusDV23-GRCz11"
ALPHA2_MINUS = "syne1a-alpha2-pred-minusDV23-GRCz11"
BIOPROJECT = "PRJNA1257549"
ENA_FASTQ = "https://ftp.sra.ebi.ac.uk/vol1/fastq"
K_CAPTURE = 20
FIRST_CAPTURED_EXON = 126          # giant exon numbering (ENSDART00000139350.4)
MIN_ANCHOR = 8                     # nt aligned on each side of a counted junction
MIN_ANCHOR_NOVEL = 12
MIN_ALIGNED_FRACTION = 0.85        # of the read length
MAX_NM_RATE = 0.05                 # edit distance / aligned length
RNADIR = DATA / "rnaseq"
HITS = RNADIR / "hits"
PATTERNS = RNADIR / "capture_kmers.txt"
PARTS = DATA / "cache" / "fastq_parts"
OUT = RESULTS / "rnaseq"

# Wild-type libraries of PRJNA1257549 (sample names from SRA RunInfo; checked in `scan`).
# run, SRA sample name, tissue, fish id. run1/run2 of one fish are the same library sequenced twice.
SAMPLES = [
    ("SRR33394562", "fast_mb_wt_42_run1", "fast muscle", "42"),
    ("SRR33394549", "slow_mb_wt_42_run1", "slow muscle", "42"),
    ("SRR33394560", "fast_mb_wt_43", "fast muscle", "43"),
    ("SRR33394546", "slow_mb_wt_43", "slow muscle", "43"),
    ("SRR33394558", "fast_mb_wt_46_run1", "fast muscle", "46"),
    ("SRR33394545", "slow_mb_wt_46", "slow muscle", "46"),
    ("SRR33394556", "fast_mb_wt_48", "fast muscle", "48"),
    ("SRR33394561", "fast_mb_wt_42_run2", "fast muscle", "42"),
    ("SRR33394557", "fast_mb_wt_46_run2", "fast muscle", "46"),
    ("SRR33394547", "slow_mb_wt_42_run2", "slow muscle", "42"),
    ("SRR33394541", "ventricle_mb_wt_5", "ventricle", "5"),
    ("SRR33394540", "ventricle_mb_wt_6", "ventricle", "6"),
    ("SRR33394539", "ventricle_mb_wt_7", "ventricle", "7"),
    ("SRR33394569", "ventricle_mb_wt_8", "ventricle", "8"),
]
TISSUES = ["fast muscle", "slow muscle", "ventricle"]
MOTIF_RANK = {"GT-AG": 0, "GC-AG": 1, "AT-AC": 2}


# ------------------------------------------------------------------------------------------------
# Locus, gene model and features
# ------------------------------------------------------------------------------------------------
def load_locus():
    built = load_json(DATA / "models_built.json")[VIEW]
    region = built["region"]
    seq = read_fasta((DATA / "genomic" / f"{VIEW}.fa").read_text())
    seq = next(iter(seq.values())).upper()
    assert len(seq) == region["end"] - region["start"] + 1, "locus FASTA does not match region"
    tx = {t["id"]: t for t in built["transcripts"]}
    giant = [tuple(e) for e in tx[GIANT]["exons"]]
    plus = [tuple(e) for e in tx[ALPHA2_PLUS]["exons"]]
    minus = [tuple(e) for e in tx[ALPHA2_MINUS]["exons"]]
    alpha2 = plus[0]
    assert alpha2 not in giant
    dv23 = [e for e in plus if e not in minus]
    assert len(dv23) == 1
    dv23 = dv23[0]
    assert plus[1] in giant, "exon after 1alpha2 must be a giant exon"
    return dict(region=region, seq=seq, start=region["start"], giant=giant, alpha2=alpha2, dv23=dv23,
                alpha2_cds_start=tx[ALPHA2_PLUS]["cds"][0], alpha2_exons=plus,
                other_tx={k: [tuple(e) for e in v["exons"]] for k, v in tx.items()})


def features(L):
    """Named exons and key junctions (1-based genomic coordinates, + strand)."""
    giant = L["giant"]
    num = {e: i + 1 for i, e in enumerate(giant)}
    i128 = giant.index(L["alpha2_exons"][1])           # exon that follows 1alpha2 (giant exon 128)
    e_after = giant[i128]
    e_before = giant[i128 - 1]                          # giant exon 127
    idv = giant.index(L["dv23"])
    a2 = L["alpha2"]
    key = {
        "alpha2": dict(label=f"1\u03b12>e{num[e_after]}", donor=a2[1], acceptor=e_after[0]),
        "giant": dict(label=f"e{num[e_before]}>e{num[e_after]}", donor=e_before[1], acceptor=e_after[0]),
        "dv23_in1": dict(label=f"e{idv}>e{idv + 1}(DV23)", donor=giant[idv - 1][1], acceptor=L["dv23"][0]),
        "dv23_in2": dict(label=f"e{idv + 1}(DV23)>e{idv + 2}", donor=L["dv23"][1], acceptor=giant[idv + 1][0]),
        "dv23_skip": dict(label=f"e{idv}>e{idv + 2} (DV23 skipped)", donor=giant[idv - 1][1],
                          acceptor=giant[idv + 1][0]),
    }
    # every annotated junction in the locus (any transcript) gets a label
    labels = {}
    for i in range(len(giant) - 1):
        labels[(giant[i][1], giant[i + 1][0])] = f"e{i + 1}>e{i + 2}"
    for k, v in key.items():
        labels[(v["donor"], v["acceptor"])] = v["label"]
    annotated = set(labels)
    for ex in L["other_tx"].values():
        for i in range(len(ex) - 1):
            annotated.add((ex[i][1], ex[i + 1][0]))
    return dict(num=num, e127=e_before, e128=e_after, n127=num[e_before], n128=num[e_after], key=key,
                labels=labels, annotated=annotated, idv=idv + 1)


# ------------------------------------------------------------------------------------------------
# patterns
# ------------------------------------------------------------------------------------------------
def _entropy(s: str) -> float:
    c = Counter(s)
    return -sum(v / len(s) * math.log2(v / len(s)) for v in c.values())


def make_patterns(L, F) -> dict:
    start, seq = L["start"], L["seq"]
    giant = L["giant"]
    captured = [e for e in giant if F["num"][e] >= FIRST_CAPTURED_EXON]
    intron127 = (F["e127"][1] + 1, F["e128"][0] - 1)
    intervals = captured + [intron127]
    lo, hi = min(a for a, _ in intervals), max(b for _, b in intervals)
    rm = ucsc_track("danRer11", "rmsk", "chr20", lo - 1, hi)
    masked = [(r["genoStart"] + 1, r["genoEnd"]) for r in rm]          # rmsk is 0-based half-open
    exonic = set()
    for a, b in captured:
        exonic.update(range(a, b + 1))

    def is_masked(a, b):  # k-mer [a, b]; exonic bases are never masked
        return any(not (b < m0 or a > m1) for m0, m1 in masked) and not all(x in exonic for x in range(a, b + 1))

    kmers, n_drop_rep, n_drop_lc = set(), 0, 0
    for a, b in intervals:
        for g in range(a, b - K_CAPTURE + 2):
            k = seq[g - start: g - start + K_CAPTURE]
            if "N" in k:
                continue
            if is_masked(g, g + K_CAPTURE - 1):
                n_drop_rep += 1
                continue
            if _entropy(k) < 1.5 or re.search(r"A{10}|C{10}|G{10}|T{10}", k):
                n_drop_lc += 1
                continue
            kmers.add(k)
            kmers.add(revcomp(k))
    RNADIR.mkdir(parents=True, exist_ok=True)
    PATTERNS.write_text("\n".join(sorted(kmers)) + "\n")
    info = dict(k=K_CAPTURE, n_patterns=len(kmers), intervals=intervals, repeats_masked=masked,
                dropped_repeat=n_drop_rep, dropped_low_complexity=n_drop_lc,
                first_exon=FIRST_CAPTURED_EXON, intron127=intron127)
    save_json(RNADIR / "capture_kmers.json", info)
    print(f"capture k-mers: {len(kmers)} (both strands) over {len(intervals)} intervals; "
          f"dropped {n_drop_rep} repeat, {n_drop_lc} low-complexity")
    return info


# ------------------------------------------------------------------------------------------------
# scan (download + filter)
# ------------------------------------------------------------------------------------------------
def ena_url(run: str, mate: int) -> str:
    digits = run[3:]
    sub = run[:6]
    d2 = (digits[6:].rjust(3, "0") + "/") if len(digits) > 6 else ""
    return f"{ENA_FASTQ}/{sub}/{d2}{run}/{run}_{mate}.fastq.gz"


def remote_size(url: str) -> int:
    for attempt in range(6):
        r = subprocess.run(["curl", "-s", "-r", "0-0", "-o", "/dev/null", "--max-time", "120",
                            "-w", "%{http_code} %header{content-range}", url], capture_output=True, text=True)
        parts = r.stdout.split()
        if parts and parts[0] == "206" and "/" in parts[-1]:
            return int(parts[-1].split("/")[-1])
        time.sleep(3 * (attempt + 1))
    raise RuntimeError(f"cannot get size of {url}: {r.stdout!r} {r.stderr!r}")


def fetch_range(url: str, a: int, b: int, dest: Path, tries: int = 8) -> None:
    want = b - a + 1
    for attempt in range(tries):
        subprocess.run(["curl", "-s", "-S", "--fail", "-r", f"{a}-{b}", "-o", str(dest), "--max-time", "1800",
                        url], capture_output=True, text=True)
        if dest.exists() and dest.stat().st_size == want:
            return
        time.sleep(4 * (attempt + 1))
    raise RuntimeError(f"range {a}-{b} of {url} failed after {tries} tries")


def runinfo() -> dict:
    ids = http(f"{EUTILS}/esearch.fcgi", params=dict(db="sra", term=f"{BIOPROJECT}[BioProject]", retmax=500,
                                                     retmode="json"), kind="ncbi", as_json=True)
    idlist = ids["esearchresult"]["idlist"]
    csv = http(f"{EUTILS}/efetch.fcgi", params=dict(db="sra", id=",".join(sorted(idlist)), rettype="runinfo",
                                                    retmode="text"), kind="ncbi")
    lines = [ln for ln in csv.splitlines() if ln.strip()]
    head = lines[0].split(",")
    rows = {}
    for ln in lines[1:]:
        if ln.startswith("Run,"):
            continue
        v = dict(zip(head, ln.split(",")))
        rows[v["Run"]] = v
    return rows


def process_file(run: str, mate: int, parts: list[Path], complete: bool, size: int, nbytes: int) -> dict:
    out = HITS / f"{run}_{mate}.tsv.gz"
    nfile = PARTS / f"{run}_{mate}.nlines"
    awk = shutil.which("mawk") or "awk"
    gz = "gzip -dc" if complete else "{ gzip -dc 2>/dev/null || true; }"
    # ripgrep's Aho-Corasick is ~3x faster than GNU grep -F with ~11,000 patterns; same output
    matcher = "rg -F --no-line-number --no-filename" if shutil.which("rg") else "LC_ALL=C grep -F"
    cmd = (f"set -o pipefail; cat {' '.join(str(p) for p in parts)} | {gz} | "
           f"{awk} -v OUT={nfile} 'NR%4==1{{h=$1}} NR%4==2{{print h \"\\t\" $0}} END{{print NR > OUT}}' | "
           f"{{ {matcher} -f {PATTERNS} || test $? -eq 1; }} | gzip -c > {out}.tmp")
    t0 = time.time()
    r = subprocess.run(["bash", "-c", cmd], capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"{run}_{mate}: filter pipeline failed: {r.stderr[-500:]}")
    nlines = int(nfile.read_text().strip())
    Path(f"{out}.tmp").rename(out)
    n_hits = sum(1 for _ in gzip.open(out, "rt"))
    info = dict(run=run, mate=mate, url=ena_url(run, mate), file_bytes=size, bytes_read=nbytes, complete=complete,
                reads_scanned=nlines // 4, reads_captured=n_hits, filter_seconds=round(time.time() - t0, 1),
                date=time.strftime("%Y-%m-%d"))
    save_json(HITS / f"{run}_{mate}.json", info)
    for p in parts:
        p.unlink(missing_ok=True)
    nfile.unlink(missing_ok=True)
    print(f"  {run}_{mate}: {info['reads_scanned']:,} reads scanned, {n_hits:,} captured "
          f"({'complete' if complete else 'PREFIX ONLY'})", flush=True)
    return info


def scan(runs: list[str], max_bytes: int | None, workers: int, chunk: int, mates=(1, 2)) -> None:
    if not PATTERNS.exists():
        raise SystemExit("run `patterns` first")
    HITS.mkdir(parents=True, exist_ok=True)
    meta = runinfo()
    known = {s[0]: s for s in SAMPLES}
    for run in runs:
        if run in known and meta.get(run, {}).get("SampleName") not in (None, known[run][1]):
            print(f"warning: {run} SampleName {meta[run].get('SampleName')} != {known[run][1]}")
    todo = []
    for run in runs:
        for mate in mates:
            info_path = HITS / f"{run}_{mate}.json"
            if info_path.exists():
                old = load_json(info_path)
                if old.get("complete") or (max_bytes and old.get("bytes_read", 0) >= max_bytes):
                    print(f"  {run}_{mate}: already scanned, skipping")
                    continue
            todo.append((run, mate))
    with ThreadPoolExecutor(workers) as dl, ThreadPoolExecutor(2) as proc:
        pending = []
        for run, mate in todo:
            url = ena_url(run, mate)
            size = remote_size(url)
            nbytes = min(size, max_bytes) if max_bytes else size
            d = PARTS / f"{run}_{mate}"
            d.mkdir(parents=True, exist_ok=True)
            ranges = [(a, min(a + chunk, nbytes) - 1) for a in range(0, nbytes, chunk)]
            parts = [d / f"part_{i:05d}" for i in range(len(ranges))]
            t0 = time.time()
            futs = [dl.submit(fetch_range, url, a, b, p) for (a, b), p in zip(ranges, parts)]
            for f in futs:
                f.result()
            dt = time.time() - t0
            print(f"{run}_{mate}: downloaded {nbytes / 1e6:,.0f} MB in {dt / 60:.1f} min "
                  f"({nbytes / 1e6 / max(dt, 1e-9):.1f} MB/s)", flush=True)
            pending.append(proc.submit(process_file, run, mate, parts, nbytes == size, size, nbytes))
        for p in pending:
            p.result()


# ------------------------------------------------------------------------------------------------
# analyse
# ------------------------------------------------------------------------------------------------
def shift_junction(seq: str, d: int, a: int) -> tuple[int, int, str]:
    """d = 0-based index of the last exon base, a = index of the first base of the next exon.
    Return the equivalent position with the best splice motif (GT-AG > GC-AG > AT-AC > other)."""
    cands = [(d, a)]
    x, y = d, a
    while x > 0 and seq[x] == seq[y - 1]:            # slide left
        x, y = x - 1, y - 1
        cands.append((x, y))
    x, y = d, a
    while y < len(seq) - 1 and seq[x + 1] == seq[y]:  # slide right
        x, y = x + 1, y + 1
        cands.append((x, y))

    def motif(c):
        return seq[c[0] + 1: c[0] + 3] + "-" + seq[c[1] - 2: c[1]]
    best = min(cands, key=lambda c: (MOTIF_RANK.get(motif(c), 3), c[0]))
    return best[0], best[1], motif(best)


def align_read(al, ref: str, read: str):
    """Primary minimap2 alignment of one read to the locus. Junctions are shifted to the preferred
    motif among equivalent positions; each is returned as
    (donor_idx, acceptor_idx, motif, left_anchor, right_anchor, mismatches_within_12nt, query_offset)."""
    hits = [h for h in al.map(read) if h.is_primary]
    if not hits:
        return None
    h = hits[0]
    q = read if h.strand == 1 else revcomp(read)
    qpos = h.q_st if h.strand == 1 else len(read) - h.q_en
    q0, rpos = qpos, h.r_st
    blocks, mm, raw, seg_m = [], set(), [], [0]
    for length, op in h.cigar:
        if op in (0, 7, 8):                       # M / = / X
            for i in range(length):
                if q[qpos + i] != ref[rpos + i]:
                    mm.add(rpos + i)
            blocks.append((rpos, rpos + length))
            qpos += length
            rpos += length
            seg_m[-1] += length
        elif op == 1:                             # insertion in the read
            qpos += length
        elif op == 2:                             # deletion in the read
            rpos += length
        elif op == 3:                             # N = intron
            raw.append((rpos - 1, rpos + length, qpos))
            rpos += length
            seg_m.append(0)
    juncs = []
    for j, (d, a, qj) in enumerate(raw):
        d2, a2, motif = shift_junction(ref, d, a)
        sh = d2 - d
        near = sum(1 for x in range(d2 - 11, d2 + 1) if x in mm) + sum(1 for x in range(a2, a2 + 12) if x in mm)
        juncs.append((d2, a2, motif, seg_m[j] + sh, seg_m[j + 1] - sh, near, qj + sh))
    aligned_q = h.q_en - h.q_st
    ok = aligned_q >= MIN_ALIGNED_FRACTION * len(read) and h.NM <= MAX_NM_RATE * aligned_q
    return dict(strand=h.strand, q=q, q0=q0, q1=qpos, blocks=blocks, juncs=juncs, NM=h.NM, ok=ok)


def read_hits(run: str):
    for mate in (1, 2):
        p = HITS / f"{run}_{mate}.tsv.gz"
        if not p.exists():
            continue
        with gzip.open(p, "rt") as fh:
            for line in fh:
                name, seq = line.rstrip("\n").split("\t")
                frag = name.lstrip("@").split()[0]
                if frag.endswith("/1") or frag.endswith("/2"):
                    frag = frag[:-2]
                yield mate, frag, seq.upper()


COUNT_FIELDS = ["reads_scanned", "reads_captured", "reads_mapped", "reads_mapped_duplicate",
                "alpha2_junction", "alpha2_junction_raw", "alpha2_junction_24mer", "alpha2_junction_24mer_raw",
                "giant_junction", "giant_junction_raw", "giant_junction_24mer", "giant_junction_24mer_raw",
                "dv23_in1", "dv23_in2", "dv23_skip", "alpha2_rescued", "giant_rescued",
                "junctions_into_alpha2_region", "unspliced_through_alpha2_donor", "unspliced_through_alpha2_5p_end"]


def rescue_clipped(r, seq: str, cands: dict, acc: int, min_clip: int = 8, slack: int = 5, ctx: int = 20):
    """Symmetric rescue of short overhangs at the exon-128 acceptor that the aligner soft-clipped instead of
    splicing (minimap2 scores non-canonical GC-AG introns lower than GT-AG ones).
    cands = {key: donor_idx} for the junctions that enter the same acceptor `acc` (0-based locus indices).
    Case 1: the alignment starts within `slack` nt of the acceptor and >= min_clip read bases are unaligned
    before it - the clipped bases are compared with the end of each candidate donor exon and with the
    unspliced intron. Case 2: the alignment ends within `slack` nt of a candidate donor and >= min_clip bases
    are unaligned after it - they are compared with the start of exon 128 and with the unspliced intron.
    The clipped bases plus `ctx` aligned bases are scored; a junction is assigned only with <= 1 mismatch and
    strictly fewer mismatches than every alternative. Returns the key or None."""
    if r["juncs"]:
        return None
    q, q0, q1 = r["q"], r["q0"], r["q1"]
    b0, last = r["blocks"][0][0], r["blocks"][-1][1] - 1
    n = len(seq)
    scores = []

    def count(idx, rng):
        m = 0
        for i in rng:
            g = idx(i)
            if g < 0 or g >= n:
                return 99
            m += q[i] != seq[g]
        return m
    if q0 >= min_clip and acc - slack <= b0 <= acc + slack:
        rng = range(0, min(len(q), q0 + ctx))
        scores.append((count(lambda i: b0 + i - q0, rng), "_unspliced"))
        for k, d in cands.items():
            scores.append((count(lambda i, d=d: (b0 + i - q0) if b0 + i - q0 >= acc else d - (acc - 1 - (b0 + i - q0)),
                                 rng), k))
    elif len(q) - q1 >= min_clip:
        rng = range(max(0, q1 - ctx), len(q))
        for k, d in cands.items():
            if d - slack <= last <= d + slack:
                scores.append((count(lambda i: last + i - (q1 - 1), rng), "_unspliced"))
                scores.append((count(lambda i, d=d: (last + i - (q1 - 1)) if last + i - (q1 - 1) <= d
                                     else acc + (last + i - (q1 - 1)) - d - 1, rng), k))
    if not scores:
        return None
    scores.sort()
    best = scores[0]
    if best[0] <= 1 and best[1] != "_unspliced" and all(x[0] > best[0] for x in scores[1:]):
        return best[1]
    return None


def analyse() -> None:
    import mappy  # minimap2 Python binding
    import numpy as np

    L = load_locus()
    F = features(L)
    start, seq = L["start"], L["seq"]
    OUT.mkdir(parents=True, exist_ok=True)
    al = mappy.Aligner(seq=seq, preset="splice:sr")
    if not al:
        raise SystemExit("minimap2 index failed")

    def norm(d, a):  # genomic donor/acceptor -> normalised genomic donor/acceptor + motif
        d2, a2, m = shift_junction(seq, d - start, a - start)
        return d2 + start, a2 + start, m
    key = {}
    for k, v in F["key"].items():
        d, a, m = norm(v["donor"], v["acceptor"])
        key[k] = dict(v, donor_norm=d, acceptor_norm=a, motif=m)
    labels = {norm(d, a)[:2]: lab for (d, a), lab in F["labels"].items()}
    annotated = {norm(d, a)[:2] for d, a in F["annotated"]}
    j24 = {}
    for k, v in key.items():
        d, a = v["donor_norm"] - start, v["acceptor_norm"] - start
        s = seq[d - 11: d + 1] + seq[a: a + 12]
        j24[k] = (s, revcomp(s))
    a2s, a2e = L["alpha2"][0] - start, L["alpha2"][1] - start
    win0, win1 = F["e127"][0] - start - 300, F["e128"][1] - start + 300      # coverage window (index)
    exon_len = {F["num"][e]: (e[0] - start, e[1] - start) for e in L["giant"] if F["num"][e] >= FIRST_CAPTURED_EXON}
    exon_len["1a2"] = (a2s, a2e)
    a2key = (key["alpha2"]["donor_norm"], key["alpha2"]["acceptor_norm"])
    resc_cands = {"alpha2": key["alpha2"]["donor_norm"] - start, "giant": key["giant"]["donor_norm"] - start}
    resc_acc = key["alpha2"]["acceptor_norm"] - start
    assert resc_acc == key["giant"]["acceptor_norm"] - start

    samples, junc_counts, cov, examples = [], defaultdict(Counter), {}, []
    present = [s for s in SAMPLES if any((HITS / f"{s[0]}_{m}.json").exists() for m in (1, 2))]
    if not present:
        raise SystemExit("no scanned runs found; run `scan` first")
    for run, name, tissue, fish in present:
        t0 = time.time()
        infos = [load_json(HITS / f"{run}_{m}.json") for m in (1, 2) if (HITS / f"{run}_{m}.json").exists()]
        diff = np.zeros(win1 - win0 + 2, dtype=np.int64)
        exon_bases = Counter()
        frag_j, frag_24 = defaultdict(set), defaultdict(set)
        frag_ju, frag_24u = defaultdict(set), defaultdict(set)   # duplicates removed
        rescued = defaultdict(set)
        seen_aln, seen_seq = set(), set()
        n_cap = n_map = n_dup = 0
        strand_by_mate = Counter()
        into_a2, retained_donor, through_5p = set(), set(), set()
        for mate, frag, rd in read_hits(run):
            n_cap += 1
            sig_seq = (mate, min(rd, revcomp(rd)))
            dup_seq = sig_seq in seen_seq
            seen_seq.add(sig_seq)
            for k, (s1, s2) in j24.items():
                if s1 in rd or s2 in rd:
                    frag_24[k].add(frag)
                    if not dup_seq:
                        frag_24u[k].add(frag)
            r = align_read(al, seq, rd)
            if not r or not r["ok"]:
                continue
            n_map += 1
            # PCR/optical duplicate: same mate, strand, alignment ends and junctions as an earlier read
            sig = (mate, r["strand"], r["blocks"][0][0], r["blocks"][-1][1], tuple(j[:2] for j in r["juncs"]))
            dup = sig in seen_aln
            seen_aln.add(sig)
            n_dup += dup
            strand_by_mate[(mate, r["strand"])] += 1
            for b0, b1 in r["blocks"]:                       # b1 exclusive
                c0, c1 = max(b0, win0), min(b1, win1 + 1)
                if c0 < c1:
                    diff[c0 - win0] += 1
                    diff[c1 - win0] -= 1
                for num, (x0, x1) in exon_len.items():
                    o = min(b1, x1 + 1) - max(b0, x0)
                    if o > 0:
                        exon_bases[num] += o
                if b0 <= a2e - 7 and b1 >= a2e + 1 + 8:         # contiguous across the 1alpha2 donor
                    retained_donor.add(frag)
                if b0 <= a2s - 8 and b1 >= a2s + 8:             # contiguous across the predicted 5' end
                    through_5p.add(frag)
            for d, a, motif, left, right, near, qj in r["juncs"]:
                g = (d + start, a + start)
                need = MIN_ANCHOR if g in annotated or g == a2key else MIN_ANCHOR_NOVEL
                if min(left, right) < need or near > 1:
                    continue
                frag_j[(g[0], g[1], motif)].add(frag)
                if not dup:
                    frag_ju[(g[0], g[1], motif)].add(frag)
                if a2s - 600 <= a <= a2e + 1:
                    into_a2.add((g, frag))
                if g == a2key and len(examples) < 400:
                    examples.append((name, mate, frag, r, qj))
            rk = rescue_clipped(r, seq, resc_cands, resc_acc)
            if rk:
                v = key[rk]
                jk = (v["donor_norm"], v["acceptor_norm"], v["motif"])
                frag_j[jk].add(frag)
                rescued[rk].add(frag)
                if not dup:
                    frag_ju[jk].add(frag)
        for jk, fr in frag_j.items():
            junc_counts[jk][run] = len(frag_ju.get(jk, ()))
            junc_counts[jk][run + "_raw"] = len(fr)
        cov[run] = np.cumsum(diff)[:win1 - win0 + 1]

        def jc(k, src=frag_j):
            v = key[k]
            return len(src.get((v["donor_norm"], v["acceptor_norm"], v["motif"]), ()))
        row = dict(run=run, sample=name, tissue=tissue, fish=fish, mates=len(infos),
                   complete=len(infos) == 2 and all(i["complete"] for i in infos),
                   reads_scanned=sum(i["reads_scanned"] for i in infos), reads_captured=n_cap, reads_mapped=n_map,
                   reads_mapped_duplicate=n_dup,
                   alpha2_junction=jc("alpha2", frag_ju), alpha2_junction_raw=jc("alpha2"),
                   alpha2_junction_24mer=len(frag_24u["alpha2"]), alpha2_junction_24mer_raw=len(frag_24["alpha2"]),
                   giant_junction=jc("giant", frag_ju), giant_junction_raw=jc("giant"),
                   giant_junction_24mer=len(frag_24u["giant"]), giant_junction_24mer_raw=len(frag_24["giant"]),
                   dv23_in1=jc("dv23_in1", frag_ju), dv23_in2=jc("dv23_in2", frag_ju),
                   dv23_skip=jc("dv23_skip", frag_ju),
                   alpha2_rescued=len(rescued["alpha2"]), giant_rescued=len(rescued["giant"]),
                   junctions_into_alpha2_region=len({x[1] for x in into_a2 if x[0] != a2key}),
                   unspliced_through_alpha2_donor=len(retained_donor),
                   unspliced_through_alpha2_5p_end=len(through_5p),
                   exon_mean_cov={str(k): round(exon_bases[k] / (x1 - x0 + 1), 2) for k, (x0, x1) in exon_len.items()},
                   strand_by_mate={f"R{m}{'+' if s == 1 else '-'}": n for (m, s), n in sorted(strand_by_mate.items())},
                   seconds=round(time.time() - t0, 1))
        samples.append(row)
        print(f"{run} {name}: captured {n_cap:,}, mapped {n_map:,}; alpha2 junction {row['alpha2_junction']} "
              f"(24-mer {row['alpha2_junction_24mer']}), giant junction {row['giant_junction']} "
              f"(24-mer {row['giant_junction_24mer']}); DV23 in {row['dv23_in1']}/{row['dv23_in2']} "
              f"skip {row['dv23_skip']}; {row['seconds']} s", flush=True)
    write_results(L, F, key, labels, annotated, samples, junc_counts, cov, (win0, win1), examples)


def per_fish(samples):
    """Pool runs of the same fish and tissue (run1 + run2 are one library)."""
    acc, runs = defaultdict(Counter), defaultdict(list)
    for s in samples:
        k = (s["tissue"], s["fish"])
        for f in COUNT_FIELDS:
            acc[k][f] += s[f]
        runs[k].append(s["run"])
    rows = []
    for (tissue, f), c in acc.items():
        a, g = c["alpha2_junction"], c["giant_junction"]
        inc = (c["dv23_in1"] + c["dv23_in2"]) / 2
        rows.append(dict(tissue=tissue, fish=f, runs=runs[(tissue, f)], **c,
                         alpha2_fraction=round(a / (a + g), 3) if a + g else None,
                         dv23_psi=round(inc / (inc + c["dv23_skip"]), 3) if inc + c["dv23_skip"] else None))
    rows.sort(key=lambda r: (TISSUES.index(r["tissue"]), int(r["fish"])))
    return rows


def wilson(k, n, z=1.96):
    if n == 0:
        return [None, None]
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return [round(max(0.0, c - h), 3), round(min(1.0, c + h), 3)]


def write_results(L, F, key, labels, annotated, samples, junc_counts, cov, win, examples):
    import numpy as np
    start, seq = L["start"], L["seq"]
    win0, win1 = win
    runs = [s["run"] for s in samples]
    cols = ["run", "sample", "tissue", "fish", "mates", "complete"] + COUNT_FIELDS
    with open(OUT / "samples.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for s in samples:
            fh.write("\t".join(str(s[c]) for c in cols) + "\n")
    fish = per_fish(samples)
    with open(OUT / "per_fish.tsv", "w") as fh:
        c2 = ["tissue", "fish", "runs"] + COUNT_FIELDS + ["alpha2_fraction", "dv23_psi"]
        fh.write("\t".join(c2) + "\n")
        for r in fish:
            fh.write("\t".join(",".join(r[c]) if isinstance(r[c], list) else str(r[c]) for c in c2) + "\n")
    tissue = []
    for t in TISSUES:
        rows = [r for r in fish if r["tissue"] == t]
        if not rows:
            continue
        a = sum(r["alpha2_junction"] for r in rows)
        g = sum(r["giant_junction"] for r in rows)
        inc = sum((r["dv23_in1"] + r["dv23_in2"]) / 2 for r in rows)
        sk = sum(r["dv23_skip"] for r in rows)
        nreads = sum(r["reads_scanned"] for r in rows)
        tissue.append(dict(
            tissue=t, fish=len(rows), fish_ids=[r["fish"] for r in rows], reads_scanned=nreads,
            alpha2_junction=a, giant_junction=g,
            alpha2_junction_24mer=sum(r["alpha2_junction_24mer"] for r in rows),
            giant_junction_24mer=sum(r["giant_junction_24mer"] for r in rows),
            alpha2_fraction_pooled=round(a / (a + g), 3) if a + g else None, alpha2_fraction_ci95=wilson(a, a + g),
            alpha2_fraction_per_fish=[r["alpha2_fraction"] for r in rows],
            alpha2_junction_per_10M_reads=round(a / nreads * 1e7, 2) if nreads else None,
            giant_junction_per_10M_reads=round(g / nreads * 1e7, 2) if nreads else None,
            dv23_psi_pooled=round(inc / (inc + sk), 3) if inc + sk else None,
            dv23_psi_per_fish=[r["dv23_psi"] for r in rows],
            junctions_into_alpha2_region=sum(r["junctions_into_alpha2_region"] for r in rows),
            unspliced_through_alpha2_donor=sum(r["unspliced_through_alpha2_donor"] for r in rows),
            unspliced_through_alpha2_5p_end=sum(r["unspliced_through_alpha2_5p_end"] for r in rows),
            exon_mean_cov={k: round(sum(s["exon_mean_cov"][k] for s in samples if s["tissue"] == t), 1)
                           for k in samples[0]["exon_mean_cov"]}))
    # ---- all junctions --------------------------------------------------------------------------
    jrows = []
    for (d, a, motif), c in junc_counts.items():
        jrows.append(dict(donor=d, acceptor=a, intron_len=a - d - 1, motif=motif, label=labels.get((d, a), ""),
                          annotated=(d, a) in annotated, total=sum(c[r] for r in runs),
                          total_raw=sum(c[r + "_raw"] for r in runs),
                          per_tissue={t: sum(c[s["run"]] for s in samples if s["tissue"] == t) for t in TISSUES},
                          per_run={r: c[r] for r in runs}))
    jrows.sort(key=lambda r: (r["donor"], r["acceptor"]))
    with open(OUT / "junctions.tsv", "w") as fh:
        fh.write("chrom\tdonor_last_exon_base\tacceptor_first_exon_base\tintron_len\tmotif\tlabel\tannotated\t"
                 "total_dedup\ttotal_raw\t"
                 + "\t".join(TISSUES) + "\t" + "\t".join(runs) + "\n")
        for r in jrows:
            fh.write(f"chr20\t{r['donor']}\t{r['acceptor']}\t{r['intron_len']}\t{r['motif']}\t{r['label']}\t"
                     f"{r['annotated']}\t{r['total']}\t{r['total_raw']}\t"
                     + "\t".join(str(r["per_tissue"][t]) for t in TISSUES)
                     + "\t" + "\t".join(str(r["per_run"][x]) for x in runs) + "\n")
    tot = sum(r["total"] for r in jrows)
    a2key = (key["alpha2"]["donor_norm"], key["alpha2"]["acceptor_norm"])
    ann_tot = sum(r["total"] for r in jrows if r["annotated"])
    consecutive = [labels[k] for k in labels if k in annotated]
    control = dict(
        junction_fragments_total=tot, annotated_fraction=round(ann_tot / tot, 4) if tot else None,
        alpha2_fraction_of_all=round(sum(r["total"] for r in jrows if (r["donor"], r["acceptor"]) == a2key) / tot, 4)
        if tot else None,
        novel_top=[dict(donor=r["donor"], acceptor=r["acceptor"], motif=r["motif"], total=r["total"])
                   for r in sorted(jrows, key=lambda r: -r["total"])
                   if not r["annotated"] and (r["donor"], r["acceptor"]) != a2key][:10],
        giant_junctions_detected=sum(1 for r in jrows if r["label"].startswith("e") and r["annotated"]),
        giant_junctions_in_capture=sum(1 for lab in consecutive if re.match(r"e(\d+)", lab)
                                       and int(re.match(r"e(\d+)", lab).group(1)) >= FIRST_CAPTURED_EXON))
    i0, i1 = F["e127"][1], F["e128"][0]
    near_a2 = [r for r in jrows if (i0 <= r["donor"] <= i1 or i0 <= r["acceptor"] <= i1) and r["total"] > 0]
    # ---- coverage and 5' end of exon 1alpha2 -----------------------------------------------------
    pos = np.arange(win0, win1 + 1) + start
    tcov = {t: sum(cov[s["run"]] for s in samples if s["tissue"] == t) for t in TISSUES
            if any(s["tissue"] == t for s in samples)}
    with open(OUT / "coverage_exon127_to_exon128.tsv", "w") as fh:
        fh.write("pos\t" + "\t".join(tcov) + "\n")
        for i, p in enumerate(pos):
            fh.write(f"{p}\t" + "\t".join(str(int(tcov[t][i])) for t in tcov) + "\n")
    a2s, a2e = L["alpha2"]
    masked = [tuple(x) for x in load_json(RNADIR / "capture_kmers.json")["repeats_masked"]]

    def ix(g):
        return g - start - win0
    bg_pos = [g for g in range(F["e127"][1] + 1, F["e128"][0])
              if not (a2s - 600 <= g <= a2e + 100) and not any(m0 - 20 <= g <= m1 + 20 for m0, m1 in masked)]
    five = {}
    for t, c in tcov.items():
        cds_med = float(np.median(c[ix(L["alpha2_cds_start"]): ix(a2e) + 1]))
        bg = float(np.median([c[ix(g)] for g in bg_pos]))
        thr = max(0.1 * cds_med, 3 * bg, 1.0)
        e5 = None
        if cds_med > 0:
            g = L["alpha2_cds_start"]
            while g > a2s - 600 and c[ix(g)] >= thr:
                g -= 1
            e5 = g + 1
        prof = {str(g): int(c[ix(g)]) for g in range(a2s - 100, L["alpha2_cds_start"] + 1, 10)}
        five[t] = dict(cds_median_cov=cds_med, intron127_background_median=bg,
                       exon_to_background=round(cds_med / max(bg, 0.5), 1), threshold=round(thr, 1),
                       estimated_5p_end=e5, predicted_5p_end=a2s, profile_every_10nt=prof)
    # ---- example reads ------------------------------------------------------------------------------
    with open(OUT / "alpha2_junction_examples.txt", "w") as fh:
        fh.write("Reads spanning the predicted exon 1alpha2 -> exon 128 junction (GRCz11 chr20, + strand).\n"
                 f"Donor (last base of exon 1alpha2) {a2key[0]:,}; acceptor (first base of exon 128) {a2key[1]:,}.\n"
                 "Each read is shown in gene orientation and split at the junction (exon 1alpha2 | exon 128);\n"
                 "the genome is shown below it, with bases that differ from the read in lower case.\n\n")
        shown = 0
        for name, mate, frag, r, qj in examples:
            if shown >= 12:
                break
            q, q0, q1 = r["q"], r["q0"], r["q1"]
            ql, qr = q[q0:qj], q[qj:q1]
            d, a = a2key[0] - start, a2key[1] - start
            gl, gr = seq[d - len(ql) + 1: d + 1], seq[a: a + len(qr)]
            if len(gl) != len(ql) or len(gr) != len(qr):
                continue
            gm = ("".join(b if b == x else b.lower() for b, x in zip(gl, ql)) + " | "
                  + "".join(b if b == x else b.lower() for b, x in zip(gr, qr)))
            fh.write(f"{name}  {frag}  R{mate}\n  read   {ql} | {qr}\n  genome {gm}\n\n")
            shown += 1
    summary = dict(
        dataset=dict(bioproject=BIOPROJECT, source="ENA FASTQ (ftp.sra.ebi.ac.uk)", layout="paired 2x150, poly(A)",
                     note="Wild-type libraries only; run1/run2 of the same fish are pooled."),
        method=dict(capture_k=K_CAPTURE, aligner="minimap2 (mappy) preset splice:sr, syne1a locus only",
                    min_anchor=MIN_ANCHOR, min_anchor_novel=MIN_ANCHOR_NOVEL,
                    min_aligned_fraction=MIN_ALIGNED_FRACTION, max_edit_rate=MAX_NM_RATE,
                    counting="fragments (the two mates of a fragment count once)"),
        coords=dict(assembly="GRCz11", chrom="chr20", alpha2_exon=list(L["alpha2"]),
                    alpha2_cds_start=L["alpha2_cds_start"], exon127=list(F["e127"]), exon128=list(F["e128"]),
                    dv23=list(L["dv23"]), n127=F["n127"], n128=F["n128"], dv23_number=F["idv"],
                    window=[int(pos[0]), int(pos[-1])]),
        key_junctions=key, samples=samples, per_fish=fish, per_tissue=tissue, five_prime=five,
        junctions_near_alpha2=near_a2, control=control,
        coverage={t: [int(x) for x in c[::10]] for t, c in tcov.items()}, coverage_step=10,
        coverage_start=int(pos[0]), examples_file="results/rnaseq/alpha2_junction_examples.txt")
    save_json(OUT / "rnaseq_summary.json", summary)
    (OUT / "syne1a_alpha2_rnaseq.svg").write_text(figure(summary, tcov, pos))
    print("\nper tissue:")
    for t in tissue:
        print(f"  {t['tissue']:12s} fish={t['fish']} alpha2={t['alpha2_junction']} (24-mer {t['alpha2_junction_24mer']}) "
              f"giant={t['giant_junction']} (24-mer {t['giant_junction_24mer']}) fraction={t['alpha2_fraction_pooled']} "
              f"CI95={t['alpha2_fraction_ci95']} per-fish={t['alpha2_fraction_per_fish']} "
              f"DV23 PSI={t['dv23_psi_pooled']} other-junctions-into-1a2={t['junctions_into_alpha2_region']} "
              f"unspliced donor={t['unspliced_through_alpha2_donor']} 5'={t['unspliced_through_alpha2_5p_end']}")
        print(f"      exon mean coverage: {t['exon_mean_cov']}")
    for t, v in five.items():
        print(f"  5' end {t}: " + str({k: v[k] for k in v if k != 'profile_every_10nt'}))
    print("control:", control)
    print("junctions touching intron 127:")
    for r in near_a2:
        print(f"  {r['donor']:,} -> {r['acceptor']:,} {r['motif']} {r['label'] or 'NOVEL'} total={r['total']} "
              f"{r['per_tissue']}")


# ------------------------------------------------------------------------------------------------
# figure (sashimi-style SVG)
# ------------------------------------------------------------------------------------------------
def figure(S, tcov, pos) -> str:
    import numpy as np
    W, ML, MR = 980, 150, 30
    track_h, gap = 92, 26
    tissues = [t for t in TISSUES if t in tcov]
    H = 70 + len(tissues) * (track_h + gap) + 120
    x0g, x1g = int(pos[0]), int(pos[-1])

    def X(g):
        return ML + (g - x0g) / (x1g - x0g) * (W - ML - MR)
    C = S["coords"]
    a2, e127, e128 = C["alpha2_exon"], C["exon127"], C["exon128"]
    pt = {t["tissue"]: t for t in S["per_tissue"]}
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" '
         'font-family="Helvetica, Arial, sans-serif" font-size="12">',
         f'<rect width="{W}" height="{H}" fill="#ffffff"/>',
         f'<text x="{ML}" y="24" font-size="15" font-weight="bold">syne1a exon 1\u03b12 in adult zebrafish RNA-seq '
         f'(GRCz11 chr20, + strand)</text>',
         f'<text x="{ML}" y="42" fill="#555">Coverage (wild-type reads pooled per tissue, {S["dataset"]["bioproject"]}) '
         'and splice junctions counted per fragment; arcs: giant (grey) and \u03b12 (red)</text>']
    y = 62
    for t in tissues:
        c = tcov[t]
        mx = max(1, float(np.percentile(c, 99.5)))
        o.append(f'<text x="{ML - 10}" y="{y + track_h / 2}" text-anchor="end" font-weight="bold">{t}</text>')
        o.append(f'<text x="{ML - 10}" y="{y + track_h / 2 + 15}" text-anchor="end" fill="#555" font-size="11">'
                 f'max {int(mx):,}</text>')
        base = y + track_h
        step = max(1, len(c) // 700)
        pts = [f"{X(int(pos[i])):.1f},{base - min(1, c[i] / mx) * (track_h - 30):.1f}" for i in range(0, len(c), step)]
        o.append(f'<polygon points="{X(int(pos[0])):.1f},{base} {" ".join(pts)} {X(int(pos[-1])):.1f},{base}" '
                 'fill="#8fb3d9" stroke="#3b6fa6" stroke-width="0.6"/>')
        o.append(f'<line x1="{ML}" x2="{W - MR}" y1="{base}" y2="{base}" stroke="#999" stroke-width="0.6"/>')
        p = pt[t]
        for (d, a, n, col, lab) in ((e127[1], e128[0], p["giant_junction"], "#777", "giant"),
                                    (a2[1], e128[0], p["alpha2_junction"], "#c0392b", "\u03b12")):
            if not n:
                continue
            xa, xb = X(d), X(a)
            hgt = 26 if lab == "giant" else 20
            o.append(f'<path d="M{xa:.1f},{y + 30} C{xa:.1f},{y + 30 - hgt} {xb:.1f},{y + 30 - hgt} {xb:.1f},{y + 30}" '
                     f'fill="none" stroke="{col}" stroke-width="{1.2 + min(4, math.log10(n + 1) * 1.5):.1f}"/>')
            o.append(f'<text x="{(xa + xb) / 2:.1f}" y="{y + 30 - hgt * 0.75 - 4:.1f}" text-anchor="middle" '
                     f'fill="{col}" font-weight="bold">{lab}: {n:,}</text>')
        y += track_h + gap
    # gene model
    gy = y + 10
    o.append(f'<line x1="{X(e127[0]):.1f}" x2="{X(e128[1]):.1f}" y1="{gy + 10}" y2="{gy + 10}" stroke="#333"/>')
    for (s, e, lab, col) in ((e127[0], e127[1], f"exon {C['n127']}", "#555"),
                             (e128[0], e128[1], f"exon {C['n128']}", "#555"),
                             (a2[0], a2[1], "exon 1\u03b12", "#c0392b")):
        o.append(f'<rect x="{X(s):.1f}" y="{gy}" width="{max(2, X(e) - X(s)):.1f}" height="20" fill="{col}"/>')
        o.append(f'<text x="{(X(s) + X(e)) / 2:.1f}" y="{gy + 36}" text-anchor="middle">{lab}</text>')
        o.append(f'<text x="{(X(s) + X(e)) / 2:.1f}" y="{gy + 50}" text-anchor="middle" fill="#555" font-size="10">'
                 f'{s:,}-{e:,}</text>')
    o.append(f'<text x="{ML}" y="{gy + 78}" fill="#555" font-size="11">Exon 1\u03b12 lies in intron {C["n127"]} of the '
             'giant transcript (ENSDART00000139350.4 numbering). Fractions and DV23 inclusion: '
             'results/rnaseq/per_fish.tsv.</text>')
    o.append("</svg>")
    return "\n".join(o)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("patterns")
    sp = sub.add_parser("scan")
    sp.add_argument("--runs", nargs="*", default=[s[0] for s in SAMPLES])
    sp.add_argument("--max-bytes", type=int, default=None, help="pilot: read only this many bytes of each file")
    sp.add_argument("--workers", type=int, default=6, help="parallel HTTP range requests")
    sp.add_argument("--chunk", type=int, default=32 * 1024 * 1024)
    sp.add_argument("--mates", type=int, nargs="*", default=[1, 2])
    sub.add_parser("analyse")
    args = ap.parse_args()
    if args.cmd == "patterns":
        L = load_locus()
        make_patterns(L, features(L))
    elif args.cmd == "scan":
        scan(args.runs, args.max_bytes, args.workers, args.chunk, tuple(args.mates))
    elif args.cmd == "analyse":
        analyse()


if __name__ == "__main__":
    main()
