#!/usr/bin/env python3
"""Step 8 - assemble the self-contained interactive viewer (SYNE1_isoform_browser.html).

Embeds the genomic loci, every transcript model, proteins, protein domains, the precomputed
alignments, the zebrafish alpha2 analysis and the c.23560G>T verification as one JSON block inside
viewer/template.html, and inlines the variant engine (viewer/variant_engine.js). The same data are
written to data/viewer_data.json for the command-line variant report (scripts/09_variant_report.js).
"""
from __future__ import annotations

import datetime as dt
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DATA, ENSEMBL, RESULTS, ROOT, http, load_json, translate  # noqa: E402

OUT = ROOT / "SYNE1_isoform_browser.html"
TEMPLATE = ROOT / "viewer" / "template.html"
ENGINE = ROOT / "viewer" / "variant_engine.js"
CHROM_ACC = {"human_SYNE1": "NC_000006.12", "mouse_Syne1": "NC_000076.7", "rat_Syne1": "NC_086019.1",
             "zebrafish_syne1a": "NC_007131.7", "zebrafish_syne1b": "NC_007128.7",
             "zebrafish_syne1a_GRCz12tu": "NC_133195.1", "zebrafish_syne1b_GRCz12tu": "NC_133192.1"}
GIANT_MSA_VIEWS = {"Human_NP_892006.3": "human_SYNE1", "Mouse_NP_001334661.1": "mouse_Syne1",
                   "Rat_NP_001406778.1": "rat_Syne1", "Zebrafish_syne1a_ENSDARP_syne1a-201": "zebrafish_syne1a",
                   "Zebrafish_syne1b_ENSDARP_syne1b-202": "zebrafish_syne1b"}


def grch37_offset(region: dict) -> int | None:
    """Constant GRCh38->GRCh37 offset if the whole human locus maps as one ungapped block."""
    mp = http(f"{ENSEMBL}/map/human/GRCh38/6:{region['start']}..{region['end']}:1/GRCh37",
              params={"content-type": "application/json"}, as_json=True).get("mappings", [])
    if len(mp) != 1:
        return None
    o, t = mp[0]["original"], mp[0]["mapped"]
    if t["strand"] != 1 or (t["end"] - t["start"]) != (o["end"] - o["start"]) or o["start"] != region["start"]:
        return None
    return t["start"] - o["start"]
DV23_RE = re.compile(r"[ST]DV[MVI]IPE[ST]PE[AV]")
ZF_A2 = "syne1a-alpha2-pred-plusDV23-GRCz11"

REFERENCES = [
    dict(text="Holt I, Duong NT, Zhang Q, et al. Specific localization of nesprin-1-α2, the short isoform of nesprin-1 "
              "with a KASH domain, in developing, fetal and regenerating muscle, using a new monoclonal antibody. "
              "BMC Cell Biol. 2016;17:26.", doi="10.1186/s12860-016-0105-9"),
    dict(text="Stroud MJ, Feng W, Zhang J, et al. Nesprin 1α2 is essential for mouse postnatal viability and nuclear "
              "positioning in skeletal muscle. J Cell Biol. 2017;216:1915-1924.", doi="10.1083/jcb.201612128"),
    dict(text="Duong NT, Morris GE, Lam LT, et al. Nesprins: tissue-specific expression of epsilon and other short "
              "isoforms. PLoS One. 2014;9:e94380.", doi="10.1371/journal.pone.0094380"),
    dict(text="Rajgor D, Mellad JA, Autore F, Zhang Q, Shanahan CM. Multiple novel nesprin-1 and nesprin-2 variants "
              "act as versatile tissue-specific intracellular scaffolds. PLoS One. 2012;7:e40098.",
         doi="10.1371/journal.pone.0040098"),
    dict(text="Holt I, Fuller HR, Lam LT, et al. Nesprin-1-alpha2 associates with kinesin at myotube outer nuclear "
              "membranes, but is restricted to neuromuscular junction nuclei in adult muscle. Sci Rep. 2019;9:14202.",
         doi="10.1038/s41598-019-50728-6"),
    dict(text="Zhou C, Li C, Zhou B, et al. Novel nesprin-1 mutations associated with dilated cardiomyopathy cause "
              "nuclear envelope disruption and defects in myogenesis. Hum Mol Genet. 2017;26:2258-2276.",
         doi="10.1093/hmg/ddx116"),
]

METHODS = [
    "Gene models: Ensembl release 116 (REST API) for human GRCh38, mouse GRCm39, rat GRCr8 and zebrafish GRCz11; "
    "NCBI RefSeq transcripts with exon and CDS coordinates from NCBI Datasets product reports; RefSeq models placed on "
    "GRCz11 from the UCSC danRer11 ncbiRefSeq track. Genomic loci (gene ± 2 kb) were downloaded from Ensembl "
    "(GRCz12tu from NCBI). mRNA, CDS and protein views are derived from the genome and the exon coordinates; "
    "deposited proteins are shown where available.",
    "Mammalian nesprin-1α2 transcripts were recognised by their protein: a unique N-terminus (MVV…DLPxCxCD) joined "
    "to the shared C-terminal region (…RVKKLKE…), as described for nesprin-1α2 by Holt et al. 2016. RefSeq names "
    "NP_001334631.1 'isoform 4 (nesprin-1 alpha 2)'. Its first exon lies in the intron before human exon 131, "
    "consistent with the mouse exon 1α2 that Stroud et al. 2017 deleted.",
    "Zebrafish: the human α2 protein was aligned (Smith-Waterman, BLOSUM62) to the syne1a and syne1b giant "
    "proteins. The residue aligned to α2 residue 31 (the codon split between exon 1α2 and exon 131) identifies "
    "the homologous exon boundary (zebrafish exon 128). The intron upstream of it was scanned for every ORF that "
    "starts with ATG, stays open, and ends at a GT or GC splice donor in the reading-frame phase needed to "
    "continue into exon 128. ORFs were ranked by Smith-Waterman similarity to the human, mouse and rat α2 "
    "peptides. Whole-locus scans with a fish α2 peptide were compared with 50 shuffled-peptide controls.",
    "The predicted α2 transcript joins exon 1α2 to syne1a exons 128-144 of Ensembl syne1a-201 (with the DV23 "
    "exon; a version without DV23 is also given). Alignments: FAMSA (pyfamsa) for proteins, a centre-star "
    "semi-global alignment for exon 1α2 nucleotides, and codon alignments derived from the protein alignment. "
    "The in-browser aligner is an affine-gap Gotoh implementation (global, local or semi-global).",
]


def load_seq(key: str) -> str:
    return "".join(open(DATA / "genomic" / f"{key}.fa").read().split("\n")[1:])


def find(view: dict, tid: str) -> dict:
    for t in view["transcripts"]:
        if t["id"] == tid or t["id"].split(".")[0] == tid.split(".")[0]:
            return t
    raise KeyError(tid)


def dv23_exons(view: dict) -> list:
    out = set()
    for t in view["transcripts"]:
        p = t.get("protein") or ""
        m = DV23_RE.search(p)
        if not m:
            continue
        r = m.start() + 3  # 1-based residue inside the exon
        for (s, e), aa in zip(t["exons"], t["exon_aa"]):
            if aa and aa[0] <= r <= aa[1] and e - s + 1 == 69:
                out.add((s, e))
    return [list(x) for x in sorted(out)]


def main() -> None:
    built = load_json(DATA / "models_built.json")
    ana = load_json(RESULTS / "alpha2_exon_analysis.json")
    alns = load_json(DATA / "alignments.json")
    primers = load_json(RESULTS / "zebrafish_alpha2_primer_suggestions.json")
    domains = load_json(DATA / "domains.json")
    vcheck_path = RESULTS / "variant_c23560G-T" / "verification.json"
    vcheck = load_json(vcheck_path) if vcheck_path.exists() else None
    sources = load_json(DATA / "sources.json")
    zx = ana["zebrafish_exon1alpha2"]
    mam = ana["mammal_exon1alpha2"]
    scans = ana["zebrafish_scans"]

    views = []
    intron_rows = []
    for key, m in built.items():
        v = dict(m["view"])
        v["region"] = m["region"]
        v["giant_ref"] = m["giant_ref"]
        v["seq"] = load_seq(key)
        assert len(v["seq"]) == m["region"]["end"] - m["region"]["start"] + 1
        giant = find(m, m["giant_ref"])
        st = giant["strand"]
        txs = []
        for t in m["transcripts"]:
            tt = {k: t.get(k) for k in ("id", "source", "name", "isoform_name", "strand", "exons", "cds", "protein",
                                        "protein_id", "protein_len", "category", "tags", "exon_aa", "predicted", "note")}
            tt = {k: val for k, val in tt.items() if val not in (None, [], "")}
            if t.get("genome_translation_matches_official") is False:
                tt["mismatch"] = (f"deposited protein differs from the genome translation at "
                                  f"{t.get('n_aa_diff_vs_genome', '?')} positions (incomplete or edited model)")
            txs.append(tt)
        # order: predicted (+DV23 first), curated RefSeq (MANE / Select first), RefSeq models, UCSC, Ensembl
        src_rank = {"Predicted (this study)": 0, "RefSeq": 1, "RefSeq (GRCz11 placement, UCSC)": 2, "Ensembl": 3}

        def tag_rank(t):
            tags = " ".join(t.get("tags", []))
            return 0 if "MANE" in tags else 1 if "REFSEQ_SELECT" in tags else 2 if "canonical" in tags else 3
        txs.sort(key=lambda t: (src_rank.get(t["source"], 4), 0 if "plusDV23" in t["id"] else 1, tag_rank(t),
                                0 if t["id"].startswith("NM_") else 1, t["id"]))
        v["transcripts"] = txs
        v["dv23"] = dv23_exons(m)
        v["domains"] = domains["views"].get(key, [])
        v["chrom_acc"] = CHROM_ACC.get(key)
        if key == "human_SYNE1":
            v["grch37_offset"] = grch37_offset(m["region"])
        # alpha2 exon highlight, hosting intron, default window
        sp = key.split("_")[0]
        if sp in mam:
            e1 = mam[sp]["exon"]
            a2 = find(m, mam[sp]["transcript"])
            k = next(i for i, ex in enumerate(giant["exons"]) if ex == a2["exons"][1])
            v["highlights"] = [dict(start=e1[0], end=e1[1], label="exon 1α2", kind="alpha2")]
            v["alpha2_status"] = "annotated"
            a2span = [min(x[0] for x in a2["exons"]), max(x[1] for x in a2["exons"])]
        else:
            fse = scans[key]["first_shared_exon"]
            k = fse["exon_number"] - 1
            if "syne1a" in key:
                e1 = zx["predicted_exon"] if key == "zebrafish_syne1a" else zx["GRCz12tu_predicted_exon"]
                v["highlights"] = [dict(start=e1[0], end=e1[1], label="exon 1α2 (predicted)", kind="alpha2")]
                v["alpha2_status"] = "predicted"
                a2 = find(m, ZF_A2 if key == "zebrafish_syne1a" else "syne1a-alpha2-pred-plusDV23-GRCz12tu")
                a2span = [min(x[0] for x in a2["exons"]), max(x[1] for x in a2["exons"])]
            else:
                ia, ib = fse["intron"]
                v["highlights"] = [dict(start=ia, end=ib, label="homologous intron, no α2 exon", kind="intron")]
                v["alpha2_status"] = "absent"
                e1 = None
                g3 = giant["exons"][-1]
                a2span = [min(ia, g3[0]), max(ib, g3[1])]
        up, down = giant["exons"][k - 1], giant["exons"][k]
        intron = [up[1] + 1, down[0] - 1] if st == 1 else [down[1] + 1, up[0] - 1]
        span = a2span[1] - a2span[0]
        v["alpha2_window"] = [max(m["region"]["start"], a2span[0] - int(span * 0.04) - 300),
                              min(m["region"]["end"], a2span[1] + int(span * 0.04) + 300)]
        v["hosting_intron"] = intron
        views.append(v)
        if key.endswith("GRCz12tu"):
            continue
        ilen = intron[1] - intron[0] + 1
        row = dict(label=f"{v['species']} {v['symbol']}", intron_len=ilen,
                   sub=f"{v['assembly']} chr{v['chrom']}:{intron[0]:,}-{intron[1]:,}",
                   up_label=f"ex {k}", down_label=f"ex {k + 1}", predicted=sp == "zebrafish")
        if e1:
            if st == 1:
                row.update(exon_from_up=e1[0] - intron[0], exon_to_down=intron[1] - e1[1])
            else:
                row.update(exon_from_up=intron[1] - e1[1], exon_to_down=e1[0] - intron[0])
            row["exon_len"] = e1[1] - e1[0] + 1
        intron_rows.append(row)

    # ---- result block -------------------------------------------------------------------
    zv = built["zebrafish_syne1a"]
    zt = find(zv, ZF_A2)
    zt_minus = find(zv, "syne1a-alpha2-pred-minusDV23-GRCz11")
    seq = load_seq("zebrafish_syne1a")
    st0 = zv["region"]["start"]
    sub = lambda a, b: seq[a - st0:b - st0 + 1]  # noqa: E731
    mrna = "".join(sub(a, b) for a, b in zt["exons"])
    cds = mrna[zt["cds"][0] - zt["exons"][0][0]:]
    cds = cds[:3 * (zt["protein_len"] + 1)]
    assert translate(cds)[:-1] == zt["protein"]
    oa = zx["other_assemblies"]
    g12 = oa["GRCz12tu (NCBI RefSeq)"]
    g13 = oa.get("GRCz13ab (NCBI RefSeq)")
    e11 = zx["predicted_exon"]
    c11 = zx["cds_part"]
    links = dict(
        ucsc=f"https://genome.ucsc.edu/cgi-bin/hgTracks?db=danRer11&position=chr20:{e11[0]}-{e11[1]}",
        ensembl=f"https://www.ensembl.org/Danio_rerio/Location/View?r=20:{e11[0]}-{e11[1]}",
        ncbi=f"https://www.ncbi.nlm.nih.gov/nuccore/NC_133195.1?report=graph&from={zx['GRCz12tu_predicted_exon'][0]}&to={zx['GRCz12tu_predicted_exon'][1]}")
    bed = (RESULTS / "zebrafish_syne1a_alpha2_predicted_GRCz11.bed").read_text()

    def fa(name, s, w=60):
        return f">{name}\n" + "\n".join(s[i:i + w] for i in range(0, len(s), w)) + "\n"
    mrna_minus = "".join(sub(a, b) for a, b in zt_minus["exons"])
    result = dict(
        grcz11=dict(exon=e11, cds=c11), coding_nt=zx["coding_nt"], peptide=zx["peptide"],
        protein_len=zt["protein_len"],
        grcz12tu=dict(acc="NC_133195.1", exon=zx["GRCz12tu_predicted_exon"], cds=[g12["start"], g12["end"]]),
        grcz13ab=dict(acc="NC_141040.1", cds=[g13["start"], g13["end"]] if g13 else None),
        links=links,
        assemblies=[
            dict(name="GRCz11 (UCSC danRer11, Ensembl 116)", acc="chr20", exon=e11, cds=c11,
                 links=[dict(label="UCSC", href=links["ucsc"]), dict(label="Ensembl", href=links["ensembl"])]),
            dict(name="GRCz12tu (current NCBI RefSeq)", acc="NC_133195.1", exon=zx["GRCz12tu_predicted_exon"],
                 cds=[g12["start"], g12["end"]], links=[dict(label="NCBI", href=links["ncbi"])]),
        ] + ([dict(name="GRCz13ab (NCBI RefSeq, AB strain)", acc="NC_141040.1", exon=None,
                   cds=[g13["start"], g13["end"]],
                   links=[dict(label="NCBI", href=f"https://www.ncbi.nlm.nih.gov/nuccore/NC_141040.1?report=graph&from={g13['start']}&to={g13['end']}")])]
             if g13 else []),
        introns=intron_rows,
        exon_utr=zx["exon_seq"][:len(zx["exon_seq"]) - zx["coding_nt"]],
        exon_cds=zx["exon_seq"][len(zx["exon_seq"]) - zx["coding_nt"]:], donor=zx["donor"],
        downloads=[
            dict(label="α2 mRNA (+DV23)", fasta=fa(f"{ZF_A2}_mRNA GRCz11 chr20 (+) {len(zt['exons'])} exons", mrna)),
            dict(label="α2 CDS", fasta=fa(f"{ZF_A2}_CDS", cds)),
            dict(label="α2 protein (+DV23)", fasta=fa(f"{ZF_A2}_protein {zt['protein_len']} aa", zt["protein"])),
            dict(label="α2 protein (−DV23)", fasta=fa(f"syne1a-alpha2-pred-minusDV23-GRCz11_protein {zt_minus['protein_len']} aa",
                                                      zt_minus["protein"])),
            dict(label="α2 mRNA (−DV23)", fasta=fa("syne1a-alpha2-pred-minusDV23-GRCz11_mRNA", mrna_minus)),
        ],
        bed=bed)

    # ---- analysis summary for the evidence tab -------------------------------------------------
    hum = [r for r in alns["alpha2_protein"]["rows"] if r[0].startswith("Human")][0][1]
    j = hum.find("TRVKKLKETL")
    end_col = j + 18 if j >= 0 else 60
    names = {"Human_NP_001334631.1": "human", "Mouse_XP_036011830.1": "mouse", "Rat_XP_063126704.1": "rat",
             "Zebrafish_syne1a_alpha2_pred": "zebrafish (pred.)", "Chicken_XP_040522257.2": "chicken",
             "Xenopus_XP_031757385.1": "Xenopus", "Herring_XP_042565831.1": "herring",
             "Pike_syne1a_XP_028970453.2": "pike syne1a", "Betta_syne1a_XP_055362111.1": "Betta syne1a",
             "Sterlet_syne1a_XP_033873712.3": "sterlet syne1a"}
    keep = ["Zebrafish_syne1a_alpha2_pred", "Herring_XP_042565831.1", "Pike_syne1a_XP_028970453.2",
            "Betta_syne1a_XP_055362111.1", "Sterlet_syne1a_XP_033873712.3", "Human_NP_001334631.1"]
    rows = {n: s for n, s in alns["alpha2_protein"]["rows"]}
    start_col = min(len(rows[n]) - len(rows[n].lstrip("-")) for n in keep)
    pep_txt = "\n".join(f"{names[n]:<17s} {rows[n][start_col:end_col]}" for n in keep)
    s1a = scans["zebrafish_syne1a"]
    s1b = scans["zebrafish_syne1b"]
    c1a = s1a["candidates"][0]
    c1b = s1b["candidates"][0] if s1b["candidates"] else None
    fse = s1a["first_shared_exon"]
    donors = " · ".join([f"{k} …{v['seq'][-6:]}|{v['donor'][:6]}" for k, v in mam.items()] +
                        [f"zebrafish …{zx['coding_seq'][-6:]}|{zx['donor'][:6]}"])
    analysis = dict(
        peptide_alignment=pep_txt,
        survey_counts=ana["vertebrate_survey_counts"], survey_total=len(ana["vertebrate_survey"]),
        scans=dict(
            syne1a=dict(first_shared=(f"syne1a exon {fse['exon_number']} ({fse['exon'][1] - fse['exon'][0] + 1} bp) = human "
                                      f"exon 131 (163 bp); split codon {fse['residue_aa']}{fse['residue']} ↔ α2 residue 31"),
                        n_candidates=s1a["n_candidates"],
                        scores=" / ".join(str(c1a["scores"][k]) for k in ("human", "mouse", "rat"))),
            syne1b=dict(intron_len=s1b["first_shared_exon"]["intron_len"],
                        scores=" / ".join(str(c1b["scores"][k]) for k in ("human", "mouse", "rat")) if c1b else "none")),
        whole_gene=dict(syne1a=ana["whole_gene_scans"]["zebrafish_syne1a"], syne1b=ana["whole_gene_scans"]["zebrafish_syne1b"]),
        donors=donors)

    # ---- alignments -------------------------------------------------------------------------------
    shorts = {"alpha2_protein": "α2 proteins, 10 vertebrates", "alpha2_cds_codon": "α2 coding sequences (codons)",
              "exon1alpha2_nt": "Exon 1α2 nucleotides", "alpha2_nterm_vertebrates": "α2 N-termini across vertebrates",
              "giant_protein": "Giant isoforms (proteins)"}
    notes = {"alpha2_protein": "first exon (unique) ends at the column before RVKKLKE",
             "exon1alpha2_nt": "annotation row: · 5'UTR, > start codon, = coding, ^ GC donor, ~ intron",
             "alpha2_nterm_vertebrates": "colour chip = group (crimson fish, amber amphibians, slate reptiles, teal birds, blue mammals)",
             "giant_protein": "8,700-8,800 residues; use the block view or search in the Sequence tab"}
    for k, s in shorts.items():
        alns[k]["short"] = s
        if k in notes:
            alns[k]["note"] = notes[k]
    # human giant exon labels for the cross-species exon table
    hv = built["human_SYNE1"]
    hg = find(hv, "NM_182961.4")
    ha2 = find(hv, "NM_001347702.2")
    for r in alns["alpha2_exon_table"]:
        ex = ha2["exons"][r["exon"] - 1]
        if r["exon"] == 1:
            r["human"]["giant"] = "unique (exon 1α2)"
        elif ex in hg["exons"]:
            r["human"]["giant"] = f"exon {hg['exons'].index(ex) + 1}"
        else:
            r["human"]["giant"] = "DV23"

    # ---- primers ------------------------------------------------------------------------------------
    primers["atg_mrna"] = zt["cds"][0] - zt["exons"][0][0] + 1
    primers["note"] = primers["note"].replace("exon 1alpha2", "exon 1\u03b12").replace("alpha2", "\u03b12")

    asm = sources["assemblies"]
    sources_text = (f"Ensembl release {sources['ensembl_release']} (" +
                    ", ".join(f"{k.replace('_', ' ')} {v['assembly_name']}" for k, v in asm.items()) +
                    "); NCBI RefSeq via the Datasets v2 API and E-utilities (zebrafish RefSeq models on GRCz12tu and "
                    "GRCz13ab); UCSC Genome Browser API (danRer11 ncbiRefSeq, mRNA and EST tracks); literature via PubMed. "
                    f"Downloaded {dt.date.today().isoformat()}.")
    data = dict(
        built_text=f"Built {dt.date.today().isoformat()} · Ensembl {sources['ensembl_release']} · NCBI RefSeq · UCSC",
        views=views, result=result, analysis=analysis, alignments=alns, primers=primers,
        methods=METHODS, sources_text=sources_text, references=REFERENCES,
        giant_msa_views=GIANT_MSA_VIEWS, domains_source=domains["source"], variant_check=vcheck)
    raw = json.dumps(data, separators=(",", ":"), ensure_ascii=False)
    (DATA / "viewer_data.json").write_text(raw)
    js = raw.replace("</", "<\\/")
    engine = ENGINE.read_text()
    assert "</script" not in engine.lower()
    html = TEMPLATE.read_text().replace("__DATA_JSON__", js).replace("/*__VARIANT_ENGINE__*/", engine)
    OUT.write_text(html)
    print(f"wrote {OUT} ({OUT.stat().st_size / 1e6:.2f} MB)")


if __name__ == "__main__":
    main()
