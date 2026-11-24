#!/usr/bin/env python3
"""Step 7 - independent verification of human SYNE1 c.23560G>T (p.Glu7854*), the homozygous
nonsense variant of the congenital muscular dystrophy patient myoblast line used by Holt et al.
2016 (BMC Cell Biol) and Gimpel et al. 2017 (Curr Biol).

Checks (each from a different source):
  1. RefSeq NM_182961.3 (used by Holt 2016 for exon numbering) and NM_182961.4 (current MANE
     Select): CDS identity, reference base at c.23560, codon 7854 and the amino-acid change.
  2. Own mapping of c.23560 onto GRCh38 through the NM_182961.4 exon model, with the genome base.
  3. Ensembl variant_recoder and VEP (RefSeq + Ensembl transcripts).
  4. GRCh38 -> GRCh37 with the Ensembl assembly mapper.
  5. Known variants overlapping codon 7854 (Ensembl) and ClinVar records (NCBI E-utilities).
Writes results/variant_c23560G-T/verification.json.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DATA, ENSEMBL, EUTILS, RESULTS, http, load_json, revcomp, save_json, translate  # noqa: E402

OUT = RESULTS / "variant_c23560G-T"
OUT.mkdir(parents=True, exist_ok=True)
C_POS, REF, ALT = 23560, "G", "T"

LITERATURE = [
    dict(ref="Holt et al. 2016, BMC Cell Biol 17:26", doi="10.1186/s12860-016-0105-9",
         says=("Immortalised myoblasts of a 16-year-old congenital muscular dystrophy patient, homozygous "
               "'SYNE1 c.23560 G>T, p.Glu7854*' (paravertebral muscle), also written 'Nucleotide: 23560 gaa to taa; "
               "Protein: E7854X'. In patient cells nesprin-1-giant mRNA was 78% lower (myoblasts and myotubes) and "
               "nesprin-1-α2 mRNA 85% lower in myotubes than the mean of five controls; neither isoform was detected "
               "by immunofluorescence; myotube formation and CKM/desmin induction were normal. The authors note that the "
               "stop codon lies 5′ to the α2 first exon and found no simple explanation for the α2 decrease."),
         check=("c.23560G>T and p.Glu7854* agree with NM_182961.3/.4. The codon is GAG, not 'gaa' as written in the "
                "methods: the change is GAG>TAG (amber stop). No known SNP makes the codon GAA, so this is most likely "
                "a typo; the protein consequence is the same. Data come from one patient line (n=1 donor) against five "
                "controls, so a line effect on α2 cannot be excluded.")),
    dict(ref="Gimpel et al. 2017, Curr Biol 27:2999", doi="10.1016/j.cub.2017.08.031",
         says=("Same immortalised patient line ('homozygous nonsense mutation within the SYNE1 gene, nucleotide "
               "23560 G>T'): nesprin-1 isoforms absent from the nuclear envelope, Pericentrin/AKAP450/PCM1 "
               "mislocalised, microtubule nucleation from the nucleus lost and nuclear spreading reduced."),
         check=("Gives no transcript version; the number only matches NM_182961 coding (c.) numbering, consistent "
                "with the verification below. The original clinical report is cited as reference 16, which the "
                "PMC full text did not include and PubMed searches for E7854X / c.23560G>T did not find.")),
    dict(ref="Attali et al. 2009, Hum Mol Genet 18:3462", doi="10.1093/hmg/ddp290",
         says="Autosomal recessive arthrogryposis caused by a SYNE1 intron 136 splice-acceptor mutation (A>G) with intron retention.",
         check="A different SYNE1 mutation; not c.23560G>T. Listed to avoid confusing the two."),
    dict(ref="Baumann et al. 2017, Eur J Hum Genet 25:262", doi="10.1038/ejhg.2016.144",
         says="Homozygous nonsense variant in the last exon of full-length SYNE1 (distal arthrogryposis); the mutant mRNA was expressed at wild-type levels.",
         check=("A different variant, but a useful contrast: a stop in the last exon escapes NMD exactly as the "
                "50-55 nt rule predicts, whereas c.23560G>T lies 16 exons upstream of the last exon.")),
]


def gb_record(acc: str) -> tuple[str, int, int, str]:
    gb = http(f"{EUTILS}/efetch.fcgi", params=dict(db="nuccore", id=acc, rettype="gb", retmode="text"), kind="ncbi")
    m = re.search(r"\n     CDS             (?:join\()?(\d+)\.\.>?(\d+)", gb)
    seq = "".join(re.findall(r"^\s+\d+ ([acgtn ]+)$", gb.split("ORIGIN")[1], re.M)).replace(" ", "").upper()
    prot = re.search(r'/translation="([^"]+)"', gb, re.S).group(1).replace("\n", "").replace(" ", "")
    return seq, int(m.group(1)), int(m.group(2)), prot


def main() -> None:
    res: dict = {"variant": "NM_182961.4:c.23560G>T", "literature": LITERATURE}
    # 1. RefSeq versions
    cds = {}
    for acc in ("NM_182961.3", "NM_182961.4"):
        seq, a, b, prot = gb_record(acc)
        c = seq[a - 1:b]
        q = (C_POS - 1) // 3
        codon = c[3 * q:3 * q + 3]
        mut = c[:C_POS - 1] + ALT + c[C_POS:]
        cds[acc] = c
        res[acc] = dict(cds=[a, b], cds_len=len(c), mrna_pos=a + C_POS - 1, ref_base=c[C_POS - 1], codon_number=q + 1,
                        codon=codon, codon_pos=(C_POS - 1) % 3 + 1, mutant_codon=mut[3 * q:3 * q + 3],
                        aa_ref=translate(codon), aa_alt=translate(mut[3 * q:3 * q + 3]), protein_aa=prot[q],
                        context=c[C_POS - 11:C_POS - 1].lower() + "[" + c[C_POS - 1] + "]" + c[C_POS:C_POS + 10].lower())
    res["cds_identical_v3_v4"] = cds["NM_182961.3"] == cds["NM_182961.4"]
    # 2. own mapping on GRCh38
    m = load_json(DATA / "models" / "human_SYNE1.json")
    t = [x for x in m["transcripts"] if x["id"] == "NM_182961.4"][0]
    g = "".join(open(DATA / "genomic" / "human_SYNE1.fa").read().split("\n")[1:])
    st = m["region"]["start"]
    start_g = t["cds"][1] if t["strand"] == -1 else t["cds"][0]
    acc_i, first, hit = 0, None, None
    for k, (s, e) in enumerate(t["exons"], 1):
        for gpos in (range(e, s - 1, -1) if t["strand"] == -1 else range(s, e + 1)):
            acc_i += 1
            if gpos == start_g:
                first = acc_i
            if first and acc_i - first + 1 == C_POS:
                hit = dict(exon=k, n_exons=len(t["exons"]), exon_coords=[s, e], exon_len=e - s + 1, grch38=gpos,
                           from_exon_5p=(e - gpos if t["strand"] == -1 else gpos - s) + 1,
                           to_exon_3p=(gpos - s if t["strand"] == -1 else e - gpos))
    plus = g[hit["grch38"] - st]
    hit.update(plus_strand_ref=plus, transcript_ref=revcomp(plus), plus_strand_alt=revcomp(ALT),
               hgvs_g=f"NC_000006.12:g.{hit['grch38']}{plus}>{revcomp(ALT)}")
    res["own_mapping"] = hit
    # 3. Ensembl variant recoder + VEP
    vr = http(f"{ENSEMBL}/variant_recoder/human/NM_182961.4:c.23560G>T", params={"content-type": "application/json"}, as_json=True)
    rec = vr[0][ALT] if vr and ALT in vr[0] else {}
    res["ensembl_variant_recoder"] = dict(hgvsg=rec.get("hgvsg"), ids=rec.get("id"),
                                          hgvsp_examples=[p for p in rec.get("hgvsp", []) if p.startswith(("NP_", "ENSP00000356224"))])
    vep = http(f"{ENSEMBL}/vep/human/hgvs/NM_182961.4:c.23560G>T",
               params={"content-type": "application/json", "refseq": 1, "numbers": 1, "hgvs": 1, "mane": 1}, as_json=True)[0]
    res["ensembl_vep"] = dict(
        most_severe=vep.get("most_severe_consequence"), seq_region=vep.get("seq_region_name"), start=vep.get("start"),
        allele_string=vep.get("allele_string"), colocated=[c.get("id") for c in vep.get("colocated_variants", [])],
        refseq_transcripts=[dict(transcript=x.get("transcript_id"), consequence=x.get("consequence_terms"), exon=x.get("exon"),
                                 codons=x.get("codons"), hgvsc=x.get("hgvsc"), hgvsp=x.get("hgvsp"))
                            for x in vep.get("transcript_consequences", []) if str(x.get("transcript_id", "")).startswith("NM_")])
    # 4. GRCh37
    mp = http(f"{ENSEMBL}/map/human/GRCh38/6:{hit['grch38']}..{hit['grch38']}:1/GRCh37", params={"content-type": "application/json"}, as_json=True)
    res["grch37"] = mp["mappings"][0]["mapped"]["start"] if mp.get("mappings") else None
    # 5. known variants in the codon + ClinVar
    p0 = hit["grch38"]
    ov = http(f"{ENSEMBL}/overlap/region/human/6:{p0 - 2}-{p0}", params={"feature": "variation", "content-type": "application/json"}, as_json=True)
    res["known_variants_in_codon"] = [dict(id=o.get("id"), pos=o.get("start"), alleles=o.get("alleles"),
                                           clinical_significance=o.get("clinical_significance"), consequence=o.get("consequence_type")) for o in ov]
    cv = {}
    for term in ['"NM_182961.4:c.23560G>T"', 'SYNE1[gene] AND "p.Glu7854Ter"', f'6[chr] AND {p0}[chrpos38]']:
        ids = http(f"{EUTILS}/esearch.fcgi", params=dict(db="clinvar", term=term, retmode="json"), kind="ncbi", as_json=True)["esearchresult"]["idlist"]
        cv[term] = ids
    codon_ids = http(f"{EUTILS}/esearch.fcgi", params=dict(db="clinvar", term="SYNE1[gene] AND E7854*", retmode="json"), kind="ncbi", as_json=True)["esearchresult"]["idlist"]
    summ = []
    if codon_ids:
        sm = http(f"{EUTILS}/esummary.fcgi", params=dict(db="clinvar", id=",".join(codon_ids), retmode="json"), kind="ncbi", as_json=True)["result"]
        for u in codon_ids:
            x = sm[u]
            cls = x.get("germline_classification") or x.get("clinical_significance") or {}
            summ.append(dict(uid=u, title=x.get("title"), classification=cls.get("description"), review=cls.get("review_status")))
    res["clinvar"] = dict(exact_searches=cv, same_codon_records=summ,
                          verdict="not in ClinVar" if not any(cv.values()) else "present in ClinVar")
    save_json(OUT / "verification.json", res)
    r4 = res["NM_182961.4"]
    print(f"NM_182961.3/.4 CDS identical: {res['cds_identical_v3_v4']}; c.23560 = {r4['ref_base']}; codon {r4['codon_number']} "
          f"{r4['codon']}>{r4['mutant_codon']} ({r4['aa_ref']}>{r4['aa_alt']})")
    print(f"GRCh38 {hit['hgvs_g']} (exon {hit['exon']}/{hit['n_exons']}); GRCh37 chr6:{res['grch37']}; "
          f"VEP {res['ensembl_vep']['most_severe']}; ClinVar: {res['clinvar']['verdict']}; same codon: {[s['title'] for s in summ]}")


if __name__ == "__main__":
    main()
