#!/usr/bin/env python3
"""Step 1 - download gene models, genomic loci and protein sequences.

Sources (all public REST APIs):
  * Ensembl REST (release reported in data/sources.json): gene/transcript/exon models,
    genomic sequence for GRCh38, GRCm39, GRCr8, GRCz11.
  * NCBI Datasets v2 product reports: RefSeq transcripts with exon coordinates and CDS.
  * NCBI E-utilities: protein FASTA, GRCz12tu genomic sequence.
  * UCSC API: ncbiRefSeq track of danRer11 (RefSeq models placed on GRCz11).

Output per view (see common.VIEWS): data/models/<key>.json and data/genomic/<key>.fa
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (DATA, FLANK, VIEWS, ens_info, ens_lookup, ens_proteins, ens_region,  # noqa: E402
                    ncbi_efetch, ncbi_product_report, ncbi_region, read_fasta, save_json,
                    ucsc_track, write_fasta)

MODELS = DATA / "models"
GENOMIC = DATA / "genomic"
for d in (MODELS, GENOMIC):
    d.mkdir(parents=True, exist_ok=True)


def tx_to_genomic(exons: list[list[int]], strand: int, tx_pos: int) -> int | None:
    """Map a 1-based transcript coordinate onto the genome through the exon list."""
    acc = 0
    for s, e in exons:
        ln = e - s + 1
        if acc + ln >= tx_pos:
            off = tx_pos - acc - 1
            return s + off if strand == 1 else e - off
        acc += ln
    return None


def ensembl_transcripts(view: dict) -> list[dict]:
    g = ens_lookup(view["ens_gene"])
    out = []
    for t in g.get("Transcript", []):
        tr = t.get("Translation")
        if t.get("biotype") != "protein_coding" or not tr:
            continue
        strand = t["strand"]
        ex = sorted(([e["start"], e["end"]] for e in t["Exon"]), key=lambda x: x[0] * strand)
        tags = []
        if t.get("is_canonical"):
            tags.append("Ensembl canonical")
        if t.get("gencode_primary"):
            tags.append("GENCODE primary")
        out.append(dict(id=f"{t['id']}.{t['version']}", source="Ensembl", name=t.get("display_name"),
                        strand=strand, exons=ex, cds=[tr["start"], tr["end"]],
                        protein_id=f"{tr['id']}.{tr['version']}", protein_len=tr["length"],
                        tags=tags, biotype=t["biotype"]))
    return out


def refseq_transcripts(view: dict, acc: str) -> list[dict]:
    rep = ncbi_product_report(view["ncbi_gene"])
    out = []
    for r in rep.get("reports", []):
        for t in r["product"].get("transcripts", []):
            loc = [l for l in t.get("genomic_locations", []) if l["genomic_accession_version"] == acc]
            if not loc:
                continue
            loc = loc[0]
            strand = 1 if loc["genomic_range"]["orientation"] == "plus" else -1
            ex = [[int(e["begin"]), int(e["end"])] for e in sorted(loc["exons"], key=lambda e: e["order"])]
            prot = t.get("protein") or {}
            cds = None
            cds_tx = None
            if t.get("cds"):
                rg = t["cds"]["range"][0]
                cds_tx = [int(rg["begin"]), int(rg["end"])]
                a = tx_to_genomic(ex, strand, cds_tx[0])
                b = tx_to_genomic(ex, strand, cds_tx[1])
                if a and b:
                    cds = sorted([a, b])
            tags = [t["select_category"]] if t.get("select_category") else []
            exon_sum = sum(e - s + 1 for s, e in ex)
            out.append(dict(id=t["accession_version"], source="RefSeq", name=t.get("name"),
                            strand=strand, exons=ex, cds=cds, cds_tx=cds_tx,
                            protein_id=prot.get("accession_version"), protein_len=prot.get("length"),
                            isoform_name=prot.get("isoform_name"), tags=tags,
                            tx_len=t.get("length"), exon_sum=exon_sum,
                            ensembl_match=t.get("ensembl_transcript"),
                            biotype=t.get("type", "").lower()))
    return out


def ucsc_refseq_transcripts(view: dict, start: int, end: int) -> list[dict]:
    items = ucsc_track(view["ucsc_db"], "ncbiRefSeq", f"chr{view['chrom']}", start, end)
    out = []
    for it in items:
        if it.get("name2") != view["symbol"]:
            continue
        strand = 1 if it["strand"] == "+" else -1
        st = [int(x) for x in str(it["exonStarts"]).strip(",").split(",")]
        en = [int(x) for x in str(it["exonEnds"]).strip(",").split(",")]
        ex = sorted(([a + 1, b] for a, b in zip(st, en)), key=lambda x: x[0] * strand)
        cds = [it["cdsStart"] + 1, it["cdsEnd"]] if it["cdsEnd"] > it["cdsStart"] else None
        out.append(dict(id=it["name"], source="RefSeq (GRCz11 placement, UCSC)", name=it["name"],
                        strand=strand, exons=ex, cds=cds, protein_id=None, protein_len=None,
                        tags=["older RefSeq annotation on GRCz11"], biotype="protein_coding"))
    return out


def main() -> None:
    sources = {"ensembl_release": ens_info("info/software")["release"],
               "assemblies": {}}
    for view in VIEWS:
        key = view["key"]
        print(f"== {key}")
        txs: list[dict] = []
        if view["ens_gene"]:
            txs += ensembl_transcripts(view)
        if view["ncbi_acc"]:
            txs += refseq_transcripts(view, view["ncbi_acc"])
        # gene span from the models collected so far
        lo = min(min(s for s, e in t["exons"]) for t in txs)
        hi = max(max(e for s, e in t["exons"]) for t in txs)
        if view.get("ucsc_refseq"):
            txs += ucsc_refseq_transcripts(view, lo - 1, hi)
            lo = min(lo, min(min(s for s, e in t["exons"]) for t in txs))
            hi = max(hi, max(max(e for s, e in t["exons"]) for t in txs))
        start, end = lo - FLANK, hi + FLANK
        if view["seq_source"] == "ensembl":
            seq = ens_region(view["ens_species"], view["chrom"], start, end, 1)
            if view["ens_species"] not in sources["assemblies"]:
                a = ens_info(f"info/assembly/{view['ens_species']}")
                sources["assemblies"][view["ens_species"]] = {
                    k: a.get(k) for k in ("assembly_name", "assembly_accession", "genebuild_last_geneset_update")}
        else:
            seq = ncbi_region(view["ncbi_acc"], start, end)
        assert len(seq) == end - start + 1, (key, len(seq), end - start + 1)
        write_fasta(GENOMIC / f"{key}.fa", [(f"{key} {view['assembly']} chr{view['chrom']}:{start}-{end}(+)", seq)])

        # proteins
        prot_ids = [t["protein_id"] for t in txs if t.get("protein_id")]
        refseq_p = [p for p in prot_ids if p[:3] in ("NP_", "XP_", "YP_")]
        ens_p = [p.split(".")[0] for p in prot_ids if p.startswith("ENS")]
        prots = {}
        if refseq_p:
            prots.update(read_fasta(ncbi_efetch("protein", refseq_p)))
        if ens_p:
            ep = ens_proteins(ens_p)
            for p in prot_ids:
                if p.startswith("ENS") and p.split(".")[0] in ep:
                    prots[p] = ep[p.split(".")[0]]
        for t in txs:
            t["protein_official"] = prots.get(t.get("protein_id")) if t.get("protein_id") else None
        model = dict(view=view, region=dict(chrom=view["chrom"], start=start, end=end), transcripts=txs)
        save_json(MODELS / f"{key}.json", model)
        n_by = {}
        for t in txs:
            n_by[t["source"]] = n_by.get(t["source"], 0) + 1
        print(f"   region chr{view['chrom']}:{start}-{end} ({end - start + 1:,} bp); transcripts {n_by}; "
              f"proteins {sum(1 for t in txs if t['protein_official'])}")
    save_json(DATA / "sources.json", sources)


if __name__ == "__main__":
    main()
