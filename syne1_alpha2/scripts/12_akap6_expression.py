#!/usr/bin/env python3
"""Step 12 - is the akap6 / nesprin-1alpha2 complex conserved across zebrafish tissues?

AKAP6 (mAKAP) is the nuclear-envelope scaffold that binds nesprin-1alpha in striated myocytes in
mammals (Pare et al. 2005, Exp Cell Res 303:388); the akap6-nesprin1alpha2 module is found in muscle
and osteoclasts (and, in mammals, heart). This step measures akap6 across zebrafish fast muscle, slow
muscle, ventricle (heart) and Ctsk+ osteoclasts, so it can be read next to where nesprin-1alpha2 (the
exon 1alpha2 -> exon 128 junction, step 10/11) is used.

Method: build gene-specific 25-mers from each gene's Ensembl cDNA, stream one+ library per tissue from
ENA, count reads assigned to a gene by a unique k-mer, and report reads/kb/million (RPKM-like). Tissue
markers verify identity (myl7/mybpc3 cardiac, mylpfa fast, smyhc1 slow, ctsk osteoclast); rpl13a
normalises; elavl3 (neuronal) is a negative control.

Result (this study): akap6 is skeletal-muscle-biased and ~14x lower in ventricle than muscle, i.e. not
expressed in the zebrafish heart - whereas nesprin-1alpha2 is cardiac-specific there. The two partners
of the mammalian complex sit in different zebrafish tissues; the cardiac mAKAP-nesprin module is not
conserved in fish.

Usage: python3 12_akap6_expression.py [--full]   (--full uses 2 libraries/tissue and no byte cap)
"""
from __future__ import annotations
import argparse, importlib.util, math, re, shutil, subprocess, sys
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from common import DATA, ENSEMBL, RESULTS, http, load_json, revcomp, save_json  # noqa: E402
_spec = importlib.util.spec_from_file_location("rj", HERE / "10_rnaseq_junctions.py")
rj = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(rj)

K = 25
CDNADIR = DATA / "cache" / "akap6_cdna"
# gene -> canonical Ensembl transcript (danio_rerio, Ensembl 116) + role
GENES = {
    "akap6":  ("ENSDART00000180856", "QUERY: nesprin-1alpha / mAKAP partner"),
    "myl7":   ("ENSDART00000040013", "cardiac marker"),
    "mybpc3": ("ENSDART00000174217", "cardiac marker"),
    "mylpfa": ("ENSDART00000075465", "fast-muscle marker"),
    "smyhc1": ("ENSDART00000171543", "slow-muscle marker"),
    "actc1a": ("ENSDART00000114549", "cardiac/skeletal actin"),
    "ctsk":   ("ENSDART00000011497", "osteoclast marker"),
    "rpl13a": ("ENSDART00000064739", "housekeeping"),
    "elavl3": ("ENSDART00000028610", "neuronal (negative control)"),
}
ORDER = ["akap6", "ctsk", "myl7", "mybpc3", "mylpfa", "smyhc1", "actc1a", "rpl13a", "elavl3"]
TISSUE_RUNS = {
    "fast muscle": ["SRR33394562", "SRR33394560"],   # PRJNA1257549 adult WT
    "slow muscle": ["SRR33394549", "SRR33394546"],
    "larval muscle (72hpf)": ["SRR26050759", "SRR26050760", "SRR26050761", "SRR26050762"],  # PRJNA1016321 mylz2+ sorted
    "ventricle":   ["SRR33394541", "SRR33394540"],
    "osteoclast":  ["SRR35522660", "SRR35522659"],   # PRJNA1331843 Ctsk+ WT sorted cells
}


def cdnas():
    CDNADIR.mkdir(parents=True, exist_ok=True)
    seqs = {}
    for g, (tx, _) in GENES.items():
        f = CDNADIR / f"{g}.fa"
        if not f.exists():
            f.write_text(http(f"{ENSEMBL}/sequence/id/{tx}", params={"type": "cdna", "content-type": "text/x-fasta"}))
        seqs[g] = "".join(l for l in f.read_text().splitlines() if not l.startswith(">")).upper()
    return seqs


def lowcomplex(k):
    if re.search(r"(A{8}|C{8}|G{8}|T{8})", k):
        return True
    c = Counter(k)
    return -sum(v / len(k) * math.log2(v / len(k)) for v in c.values()) < 1.6


def scan_run(run, mate, patt, k2g, cap):
    url = rj.ena_url(run, mate)
    size = rj.remote_size(url)
    nbytes = min(size, cap) if cap else size
    d = rj.PARTS / f"akap6_{run}_{mate}"; d.mkdir(parents=True, exist_ok=True)
    chunk = 32 * 1024 * 1024
    ranges = [(a, min(a + chunk, nbytes) - 1) for a in range(0, nbytes, chunk)]
    parts = [d / f"p{i:04d}" for i in range(len(ranges))]
    with ThreadPoolExecutor(4) as ex:
        for f in [ex.submit(rj.fetch_range, url, a, b, p) for (a, b), p in zip(ranges, parts)]:
            f.result()
    nlines = d / "nlines"
    matcher = "rg -F --no-line-number --no-filename" if shutil.which("rg") else "LC_ALL=C grep -F"
    gz = "gzip -dc" if nbytes == size else "{ gzip -dc 2>/dev/null || true; }"
    cmd = (f"set -o pipefail; cat {' '.join(str(p) for p in parts)} | {gz} | "
           f"mawk -v OUT={nlines} 'NR%4==2{{print}} END{{print NR > OUT}}' | "
           f"{{ {matcher} -f {patt} || test $? -eq 1; }}")
    r = subprocess.run(["bash", "-c", cmd], capture_output=True, text=True)
    nr = int(nlines.read_text().strip()) // 4 if nlines.exists() else 0
    for p in parts:
        p.unlink(missing_ok=True)
    nlines.unlink(missing_ok=True)
    d.rmdir()
    if r.returncode != 0:
        raise RuntimeError(r.stderr[-200:])
    per = Counter()
    for seq in r.stdout.splitlines():
        hits = {k2g[seq[i:i + K]] for i in range(len(seq) - K + 1) if seq[i:i + K] in k2g}
        if len(hits) == 1:
            per[next(iter(hits))] += 1
    return nr, per


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", action="store_true", help="2 libraries/tissue, no byte cap")
    ap.add_argument("--mates", type=int, nargs="*", default=[1])
    ap.add_argument("--tissues", nargs="*", default=None,
                    help="only (re)scan these tissues and merge into the existing table; default all")
    a = ap.parse_args()
    cap = None if a.full else 700_000_000
    runs_per = None if a.full else 1
    todo = a.tissues if a.tissues else list(TISSUE_RUNS)
    for t in todo:
        if t not in TISSUE_RUNS:
            raise SystemExit(f"unknown tissue {t!r}; choices: {list(TISSUE_RUNS)}")
    seq = cdnas()
    kmer_genes = defaultdict(set)
    for g, s in seq.items():
        for i in range(len(s) - K + 1):
            kk = s[i:i + K]
            if "N" in kk or lowcomplex(kk):
                continue
            kmer_genes[kk].add(g); kmer_genes[revcomp(kk)].add(g)
    k2g = {k: next(iter(gs)) for k, gs in kmer_genes.items() if len(gs) == 1}
    patt = CDNADIR / "panel_kmers.txt"; patt.write_text("\n".join(sorted(k2g)) + "\n")
    print("gene-specific 25-mers: " + ", ".join(f"{g}={sum(1 for v in k2g.values() if v == g)}" for g in seq))

    # merge base: keep previously-measured tissues that are not being re-scanned
    jpath = RESULTS / "zebrafish_akap6_expression.json"
    prev = load_json(jpath) if jpath.exists() else {}
    rpkm = dict(prev.get("rpkm", {})); reads = dict(prev.get("reads_scanned", {}))
    for tissue in todo:
        per = Counter(); tot = 0
        for run in TISSUE_RUNS[tissue][:runs_per]:
            for mate in a.mates:
                try:
                    nr, p = scan_run(run, mate, patt, k2g, cap)
                    tot += nr; per.update(p)
                    print(f"  [{tissue}] {run}_{mate}: {nr:,} reads, {sum(p.values()):,} on-panel", flush=True)
                except Exception as e:
                    print(f"  [{tissue}] {run}_{mate}: FAILED {str(e)[:120]}", flush=True)
        reads[tissue] = tot
        rpkm[tissue] = {g: (per[g] / (len(seq[g]) / 1000) / (tot / 1e6) if tot else 0) for g in seq}

    cols = [t for t in TISSUE_RUNS if t in rpkm]   # all tissues with data, in canonical order
    out = RESULTS / "zebrafish_akap6_expression.tsv"
    with open(out, "w") as fh:
        fh.write("gene\t" + "\t".join(cols) + "\trole\n")
        for g in ORDER:
            fh.write(f"{g}\t" + "\t".join(f"{rpkm[t].get(g, 0):.2f}" for t in cols) + f"\t{GENES[g][1]}\n")
    save_json(jpath, dict(rpkm=rpkm, reads_scanned=reads, cdna_len={g: len(s) for g, s in seq.items()},
                          genes=GENES, note="RPKM-like (reads/kb/million) from cDNA 25-mer capture; one library/tissue unless --full"))
    print(f"\n{'gene':8s} " + "".join(f"{t[:12]:>13s}" for t in cols) + "   role")
    for g in ORDER:
        print(f"{g:8s} " + "".join(f"{rpkm[t].get(g, 0):13.1f}" for t in cols) + f"   {GENES[g][1]}")
    ak = {t: rpkm[t]["akap6"] for t in cols}
    print("\nakap6 RPKM: " + ", ".join(f"{t} {v:.1f}" for t, v in ak.items()))


if __name__ == "__main__":
    main()
