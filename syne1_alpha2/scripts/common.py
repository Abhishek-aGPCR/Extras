"""Shared helpers: cached HTTP access to Ensembl / NCBI / UCSC, sequence utilities and
the gene configuration used by every step of the pipeline.

All network responses are cached under data/cache/ so that re-running a step is fast
and fully reproducible offline once the cache exists.
"""
from __future__ import annotations

import hashlib
import json
import os
import time
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
CACHE = DATA / "cache"
RESULTS = ROOT / "results"
for _d in (DATA, CACHE, RESULTS):
    _d.mkdir(parents=True, exist_ok=True)

ENSEMBL = "https://rest.ensembl.org"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
DATASETS = "https://api.ncbi.nlm.nih.gov/datasets/v2"
UCSC = "https://api.genome.ucsc.edu"

# ---------------------------------------------------------------------------
# Gene / assembly configuration
# ---------------------------------------------------------------------------
# "views" are what the HTML viewer shows: one genomic locus on one assembly.
VIEWS = [
    dict(key="human_SYNE1", species="Human", latin="Homo sapiens", symbol="SYNE1",
         assembly="GRCh38", ucsc_db="hg38", chrom="6", ncbi_acc="NC_000006.12",
         ncbi_gene=23345, ens_gene="ENSG00000131018", ens_species="homo_sapiens",
         seq_source="ensembl"),
    dict(key="mouse_Syne1", species="Mouse", latin="Mus musculus", symbol="Syne1",
         assembly="GRCm39", ucsc_db="mm39", chrom="10", ncbi_acc="NC_000076.7",
         ncbi_gene=64009, ens_gene="ENSMUSG00000096054", ens_species="mus_musculus",
         seq_source="ensembl"),
    dict(key="rat_Syne1", species="Rat", latin="Rattus norvegicus", symbol="Syne1",
         assembly="GRCr8", ucsc_db="rn8", chrom="1", ncbi_acc="NC_086019.1",
         ncbi_gene=499010, ens_gene="ENSRNOG00000018960", ens_species="rattus_norvegicus",
         seq_source="ensembl"),
    dict(key="zebrafish_syne1a", species="Zebrafish", latin="Danio rerio", symbol="syne1a",
         assembly="GRCz11", ucsc_db="danRer11", chrom="20", ncbi_acc=None,
         ncbi_gene=569553, ens_gene="ENSDARG00000009499", ens_species="danio_rerio",
         seq_source="ensembl", ucsc_refseq=True),
    dict(key="zebrafish_syne1b", species="Zebrafish", latin="Danio rerio", symbol="syne1b",
         assembly="GRCz11", ucsc_db="danRer11", chrom="17", ncbi_acc=None,
         ncbi_gene=567576, ens_gene="ENSDARG00000063068", ens_species="danio_rerio",
         seq_source="ensembl", ucsc_refseq=True),
    # Current NCBI RefSeq annotation of zebrafish is on GRCz12tu (no Ensembl models there).
    dict(key="zebrafish_syne1a_GRCz12tu", species="Zebrafish", latin="Danio rerio", symbol="syne1a",
         assembly="GRCz12tu", ucsc_db="GCF_049306965.1", chrom="20", ncbi_acc="NC_133195.1",
         ncbi_gene=569553, ens_gene=None, ens_species=None, seq_source="ncbi"),
    dict(key="zebrafish_syne1b_GRCz12tu", species="Zebrafish", latin="Danio rerio", symbol="syne1b",
         assembly="GRCz12tu", ucsc_db="GCF_049306965.1", chrom="17", ncbi_acc="NC_133192.1",
         ncbi_gene=567576, ens_gene=None, ens_species=None, seq_source="ncbi"),
]
VIEW = {v["key"]: v for v in VIEWS}
FLANK = 2000  # bp of flanking sequence kept around each gene

# ---------------------------------------------------------------------------
# Cached HTTP
# ---------------------------------------------------------------------------
_last_call: dict[str, float] = {}
_MIN_INTERVAL = {"ensembl": 0.08, "ncbi": 0.36, "ucsc": 0.3}


def _throttle(kind: str) -> None:
    now = time.time()
    wait = _MIN_INTERVAL.get(kind, 0.1) - (now - _last_call.get(kind, 0))
    if wait > 0:
        time.sleep(wait)
    _last_call[kind] = time.time()


def http(url: str, params: dict | None = None, data: dict | None = None,
         headers: dict | None = None, kind: str = "ensembl", as_json: bool = False,
         tries: int = 5):
    """GET (or POST when data is given) with an on-disk cache keyed by the request."""
    key = json.dumps([url, params, data], sort_keys=True)
    h = hashlib.sha1(key.encode()).hexdigest()
    path = CACHE / f"{h}.{'json' if as_json else 'txt'}"
    if path.exists():
        txt = path.read_text()
        return json.loads(txt) if as_json else txt
    last_err = None
    for attempt in range(tries):
        _throttle(kind)
        try:
            if data is not None:
                r = requests.post(url, params=params, data=data if not isinstance(data, str) else data,
                                  headers=headers, timeout=300)
            else:
                r = requests.get(url, params=params, headers=headers, timeout=300)
            if r.status_code == 429 or r.status_code >= 500:
                time.sleep(2 ** attempt)
                last_err = f"HTTP {r.status_code}"
                continue
            r.raise_for_status()
            txt = r.text
            if as_json:
                json.loads(txt)  # validate before caching
            path.write_text(txt)
            return json.loads(txt) if as_json else txt
        except (requests.RequestException, ValueError) as e:  # network or JSON error
            last_err = e
            time.sleep(2 ** attempt)
    raise RuntimeError(f"request failed: {url} {params} -> {last_err}")


# --- Ensembl ---------------------------------------------------------------
def ens_lookup(stable_id: str) -> dict:
    return http(f"{ENSEMBL}/lookup/id/{stable_id}",
                params={"expand": 1, "utr": 1, "content-type": "application/json"}, as_json=True)


def ens_region(species: str, chrom: str, start: int, end: int, strand: int = 1) -> str:
    """Genomic sequence (1-based inclusive). Ensembl caps a region at 10 Mb."""
    txt = http(f"{ENSEMBL}/sequence/region/{species}/{chrom}:{start}..{end}:{strand}",
               params={"content-type": "text/plain"})
    return txt.strip().upper()


def ens_proteins(translation_ids: list[str]) -> dict[str, str]:
    out = {}
    for i in range(0, len(translation_ids), 50):
        chunk = translation_ids[i:i + 50]
        res = http(f"{ENSEMBL}/sequence/id", params={"type": "protein"},
                   data=json.dumps({"ids": chunk}),
                   headers={"Content-Type": "application/json", "Accept": "application/json"},
                   as_json=True)
        for r in res:
            out[r["id"]] = r["seq"]
    return out


def ens_info(path: str) -> dict:
    return http(f"{ENSEMBL}/{path}", params={"content-type": "application/json"}, as_json=True)


# --- NCBI ------------------------------------------------------------------
def ncbi_product_report(gene_id: int) -> dict:
    return http(f"{DATASETS}/gene/id/{gene_id}/product_report", kind="ncbi", as_json=True)


def ncbi_gene_report(gene_id: int) -> dict:
    return http(f"{DATASETS}/gene/id/{gene_id}", kind="ncbi", as_json=True)


def ncbi_efetch(db: str, ids: list[str], rettype: str = "fasta", **extra) -> str:
    out = []
    for i in range(0, len(ids), 150):
        chunk = ids[i:i + 150]
        out.append(http(f"{EUTILS}/efetch.fcgi", data=dict(db=db, id=",".join(chunk), rettype=rettype,
                                                          retmode="text", **extra), kind="ncbi"))
    return "".join(out)


def ncbi_region(acc: str, start: int, end: int) -> str:
    txt = http(f"{EUTILS}/efetch.fcgi", params=dict(db="nuccore", id=acc, seq_start=start, seq_stop=end,
                                                    strand=1, rettype="fasta", retmode="text"), kind="ncbi")
    return "".join(txt.split("\n")[1:]).replace("\n", "").strip().upper()


def ncbi_esearch(db: str, term: str, retmax: int = 5000) -> list[str]:
    res = http(f"{EUTILS}/esearch.fcgi", params=dict(db=db, term=term, retmax=retmax, retmode="json"),
               kind="ncbi", as_json=True)
    return res["esearchresult"]["idlist"]


# --- UCSC ------------------------------------------------------------------
def ucsc_track(genome: str, track: str, chrom: str, start: int, end: int) -> list[dict]:
    res = http(f"{UCSC}/getData/track", params=dict(genome=genome, track=track, chrom=chrom,
                                                     start=start, end=end), kind="ucsc", as_json=True)
    items = res.get(track, [])
    if isinstance(items, dict):
        items = items.get(chrom, [])
    return items


# ---------------------------------------------------------------------------
# Sequence helpers
# ---------------------------------------------------------------------------
_COMP = str.maketrans("ACGTNacgtnRYKMSWBDHVrykmswbdhv", "TGCANtgcanYRMKSWVHDByrmkswvhdb")
CODON = {}
_b = "TCAG"
_aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for i, a in enumerate(_b):
    for j, b in enumerate(_b):
        for k, c in enumerate(_b):
            CODON[a + b + c] = _aa[16 * i + 4 * j + k]


def revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


def translate(nt: str) -> str:
    nt = nt.upper()
    return "".join(CODON.get(nt[i:i + 3], "X") for i in range(0, len(nt) - len(nt) % 3, 3))


def read_fasta(text: str) -> dict[str, str]:
    seqs, name = {}, None
    for line in text.splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif name and line.strip():
            seqs[name].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def write_fasta(path: Path, records: list[tuple[str, str]], width: int = 60) -> None:
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")


def save_json(path: Path, obj) -> None:
    path.write_text(json.dumps(obj, indent=1))


def load_json(path: Path):
    return json.loads(Path(path).read_text())
