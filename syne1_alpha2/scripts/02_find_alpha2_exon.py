#!/usr/bin/env python3
"""Step 2 - characterise the nesprin-1 alpha2 first exon (exon 1alpha2) in mammals and
locate its orthologue in zebrafish syne1a / syne1b.

Logic
-----
1. Mammalian alpha2 transcripts are recognised by their protein: a unique N-terminus
   (MVV...DLPxCxCD) joined to the shared C-terminal region (...RVKKLKE...).  Their first
   exon is exon 1alpha2; we record its coordinates, coding length, phase and splice donor.
2. In each zebrafish giant model we find the exon homologous to the first exon shared by
   alpha2 and the giant isoform in mammals (human exon 131): the human alpha2 protein is
   aligned to the zebrafish giant protein and the residue aligned to alpha2 residue 31
   (the codon split between exon 1alpha2 and exon 131) identifies the exon boundary.
3. The intron upstream of that exon is scanned for ORFs that (a) start with ATG,
   (b) stay open, (c) end at a GT or GC splice donor (d) in the reading-frame phase required to
   continue into the shared exon. Candidates are ranked by Smith-Waterman similarity to the
   mammalian alpha2 unique peptides; nucleotide similarity is also reported.
4. Controls: whole-gene peptide scans of syne1a and syne1b against shuffled peptides.
5. The exon is placed on GRCz12tu and GRCz13ab by exact sequence matching.
6. Vertebrate survey: RefSeq proteins of the alpha2 type across vertebrates.
"""
from __future__ import annotations

import random
import re
import sys
from pathlib import Path

import parasail

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (DATA, RESULTS, load_json, ncbi_efetch, ncbi_esearch, ncbi_region,  # noqa: E402
                    read_fasta, revcomp, save_json, translate, write_fasta)

MODELS = DATA / "models"
GENOMIC = DATA / "genomic"
B62 = parasail.blosum62
ALPHA2_RE = re.compile(r"C[DN]CD[A-Z]{1,2}R[VI]KKLKE")

GIANT_REF = {  # representative giant (full-length) transcript per view
    "human_SYNE1": "NM_182961.4", "mouse_Syne1": "NM_001347732.2", "rat_Syne1": "NM_001419849.1",
    "zebrafish_syne1a": "ENSDART00000139350", "zebrafish_syne1b": "ENSDART00000148520",
    "zebrafish_syne1a_GRCz12tu": "XM_073933275.1", "zebrafish_syne1b_GRCz12tu": "XM_087227815.1",
}
ALPHA2_REP = {"human_SYNE1": "NM_001347702.2", "mouse_Syne1": "NM_022027.3", "rat_Syne1": "XM_063270634.2"}


# ---------------------------------------------------------------------------
class Locus:
    def __init__(self, key: str):
        m = load_json(MODELS / f"{key}.json")
        self.key = key
        self.view = m["view"]
        self.start = m["region"]["start"]
        self.end = m["region"]["end"]
        self.seq = "".join(open(GENOMIC / f"{key}.fa").read().split("\n")[1:])
        self.tx = {t["id"]: t for t in m["transcripts"]}

    def sub(self, a: int, b: int, strand: int = 1) -> str:
        s = self.seq[a - self.start:b - self.start + 1]
        return s if strand == 1 else revcomp(s)

    def find_tx(self, prefix: str) -> dict:
        for k, t in self.tx.items():
            if k == prefix or k.split(".")[0] == prefix.split(".")[0]:
                return t
        raise KeyError(prefix)


def cds_segments(t: dict) -> list[dict]:
    """Per exon (transcript order): genomic CDS part and cumulative CDS offsets."""
    a, b = t["cds"]
    out, acc = [], 0
    for i, (s, e) in enumerate(t["exons"], 1):
        cs, ce = max(s, a), min(e, b)
        if cs <= ce:
            out.append(dict(n=i, exon=[s, e], cds=[cs, ce], start=acc, end=acc + ce - cs + 1))
            acc += ce - cs + 1
        else:
            out.append(dict(n=i, exon=[s, e], cds=None, start=acc, end=acc))
    return out


def cds_seq(loc: Locus, t: dict) -> str:
    st = t["strand"]
    return "".join(loc.sub(seg["cds"][0], seg["cds"][1], st) for seg in cds_segments(t) if seg["cds"])


def sw(q: str, r: str, go=10, ge=1, matrix=B62):
    return parasail.sw_trace_striped_16(q, r, go, ge, matrix)


def aligned_map(res) -> dict[int, int]:
    tb = res.traceback
    qi, ri, m = res.cigar.beg_query, res.cigar.beg_ref, {}
    for a, b in zip(tb.query, tb.ref):
        if a != "-" and b != "-":
            m[qi] = ri
        if a != "-":
            qi += 1
        if b != "-":
            ri += 1
    return m


# ---------------------------------------------------------------------------
def mammal_exon1(loc: Locus, tid: str) -> dict:
    t = loc.find_tx(tid)
    st = t["strand"]
    e1 = t["exons"][0]
    seq = loc.sub(e1[0], e1[1], st)
    cds_start = t["cds"][1] if st == -1 else t["cds"][0]
    utr = (e1[1] - cds_start) if st == -1 else (cds_start - e1[0])
    coding = seq[utr:]
    donor = loc.sub(e1[0] - 10, e1[0] - 1, st)[:10] if st == -1 else loc.sub(e1[1] + 1, e1[1] + 10, st)
    up = loc.sub(e1[1] + 1, e1[1] + 300, st) if st == -1 else loc.sub(e1[0] - 300, e1[0] - 1, st)
    prot = t["protein_official"]
    return dict(view=loc.key, transcript=tid, protein_id=t["protein_id"], exon=e1, strand=st,
                length=e1[1] - e1[0] + 1, utr5=utr, coding_nt=len(coding), end_phase=len(coding) % 3,
                peptide=translate(coding[:len(coding) - len(coding) % 3]), donor=donor[:8],
                seq=seq, upstream300=up, protein_nterm=prot[:40],
                second_exon=t["exons"][1])


def first_shared_exon(loc: Locus, giant: dict, human_alpha2_prot: str) -> dict:
    zp = giant["protein_official"]
    res = sw(human_alpha2_prot, zp, 11, 1)
    amap = aligned_map(res)
    # human alpha2 residue 31 (0-based 30) is the codon split between exon1alpha2 and exon 131
    zres = amap[30]
    third_nt = 3 * zres + 2  # 0-based CDS offset of the codon's third nucleotide
    segs = cds_segments(giant)
    for k, seg in enumerate(segs):
        if seg["cds"] and seg["start"] <= third_nt < seg["end"]:
            up_contrib = seg["start"] % 3  # nt of the split codon supplied by the upstream exon
            prev = segs[k - 1]
            st = giant["strand"]
            intron = ([prev["exon"][1] + 1, seg["exon"][0] - 1] if st == 1
                      else [seg["exon"][1] + 1, prev["exon"][0] - 1])
            tb = res.traceback
            ident = sum(1 for a, b in zip(tb.query, tb.ref) if a == b and a != "-")
            return dict(giant=giant["id"], exon_number=seg["n"], exon=seg["exon"], residue=zres + 1,
                        residue_aa=zp[zres], upstream_contribution=up_contrib, intron=intron,
                        intron_len=intron[1] - intron[0] + 1, prev_exon_number=prev["n"],
                        alpha2_vs_giant_identity=f"{ident}/{len(tb.query)}")
    raise RuntimeError("shared exon not found")


def scan_intron(loc: Locus, fse: dict, strand: int, peptides: dict, mam_exons: dict) -> dict:
    a, b = fse["intron"]
    intr = loc.sub(a, b, strand)
    to_g = (lambda i: a + i) if strand == 1 else (lambda i: b - i)
    need = fse["upstream_contribution"]
    cands = []
    for m in re.finditer("(?=ATG)", intr):
        s = m.start()
        for d in range(s + 30, min(len(intr) - 2, s + 900)):
            L = d - s
            if "*" in translate(intr[s:s + L - L % 3]):
                break
            if L % 3 == need and intr[d:d + 2] in ("GT", "GC"):
                pep = translate(intr[s:s + L - L % 3])
                scores = {k: sw(p, pep).score for k, p in peptides.items()}
                cands.append(dict(atg=s, donor=d, coding_nt=L, peptide=pep, donor_seq=intr[d:d + 8],
                                  donor_type=intr[d:d + 2], kozak=intr[max(0, s - 6):s + 4],
                                  scores=scores, best=max(scores.values()),
                                  genomic=sorted([to_g(s), to_g(d - 1)])))
    cands.sort(key=lambda c: (-c["best"], c["coding_nt"]))
    # nucleotide similarity of mammalian exon 1alpha2 to the intron
    nts = {}
    for k, ex in mam_exons.items():
        r = sw(ex["seq"], intr, 5, 2, parasail.dnafull)
        tb = r.traceback
        ident = sum(1 for x, y in zip(tb.query, tb.ref) if x == y and x != "-")
        nts[k] = dict(score=r.score, identity=f"{ident}/{len(tb.query)}",
                      intron_range=[r.cigar.beg_ref, r.end_ref])
    return dict(intron_seq=intr, candidates=cands[:15], n_candidates=len(cands), nt_similarity=nts)


def whole_gene_scan(loc: Locus, strand: int, pep: str, n_shuffle: int = 50, seed: int = 1) -> dict:
    s = loc.seq if strand == 1 else revcomp(loc.seq)
    frames = [translate(s[f:]) for f in range(3)]

    def best(p):
        bb = (0, None, None)
        for f, tr in enumerate(frames):
            r = sw(p, tr)
            if r.score > bb[0]:
                nt = f + 3 * r.end_ref
                g = loc.start + nt if strand == 1 else loc.end - nt
                bb = (r.score, g, r.traceback.ref)
        return bb

    sc, g, ref = best(pep)
    rnd = random.Random(seed)
    ctrl = []
    for _ in range(n_shuffle):
        p = list(pep)
        rnd.shuffle(p)
        ctrl.append(best("".join(p))[0])
    ctrl.sort()
    return dict(query=pep, best_score=sc, best_hit_genomic=g, best_hit_seq=ref,
                shuffled_mean=round(sum(ctrl) / len(ctrl), 1), shuffled_max=ctrl[-1],
                shuffled_p95=ctrl[int(0.95 * len(ctrl)) - 1])


def vertebrate_survey() -> list[dict]:
    term = ("(nesprin-1[Protein Name] OR syne1[Gene Name] OR syne1a[Gene Name] OR syne1b[Gene Name]) "
            "AND txid7742[Organism:exp] AND 700:1400[SLEN] AND refseq[filter]")
    ids = ncbi_esearch("protein", term)
    fa = read_fasta(ncbi_efetch("protein", ids))
    gp = ncbi_efetch("protein", ids, rettype="gp")
    meta = {}
    for rec in gp.split("\n//"):
        acc = re.search(r"VERSION\s+(\S+)", rec)
        if not acc:
            continue
        org = re.search(r"  ORGANISM\s+(.+?)\n((?:\s{12}.+\n)+)", rec)
        gene = re.search(r'/gene="([^"]+)"', rec)
        lineage = re.sub(r"\s+", " ", org.group(2)).strip() if org else ""
        meta[acc.group(1)] = dict(organism=org.group(1).strip() if org else "", lineage=lineage,
                                  gene=gene.group(1) if gene else "")
    groups = [("Mammalia", "Mammals"), ("Aves", "Birds"), ("Crocodylia", "Reptiles"), ("Testudines", "Reptiles"),
              ("Lepidosauria", "Reptiles"), ("Amphibia", "Amphibians"), ("Actinopterygii", "Ray-finned fish"),
              ("Coelacanthiformes", "Lobe-finned fish"), ("Dipnomorpha", "Lungfish"),
              ("Chondrichthyes", "Cartilaginous fish")]
    out = []
    for acc, s in fa.items():
        m = ALPHA2_RE.search(s[:90])
        if not m:
            continue
        md = meta.get(acc, {})
        grp = next((g for k, g in groups if k in md.get("lineage", "")), "Other vertebrates")
        uniq_end = m.start() + s[m.start():m.end()].index("R")  # residue before RVKKLKE (split codon)
        out.append(dict(accession=acc, organism=md.get("organism", ""), group=grp, gene=md.get("gene", ""),
                        length=len(s), unique_nterm=s[:uniq_end], junction=s[uniq_end:uniq_end + 12]))
    out.sort(key=lambda r: (r["group"], r["organism"]))
    return out


# ---------------------------------------------------------------------------
def main() -> None:
    loci = {k: Locus(k) for k in ["human_SYNE1", "mouse_Syne1", "rat_Syne1", "zebrafish_syne1a",
                                  "zebrafish_syne1b", "zebrafish_syne1a_GRCz12tu", "zebrafish_syne1b_GRCz12tu"]}
    report: dict = {}

    # 1. mammalian alpha2 transcripts and exon 1alpha2
    alpha2_tx = {}
    for key in ["human_SYNE1", "mouse_Syne1", "rat_Syne1"]:
        alpha2_tx[key] = [t["id"] for t in loci[key].tx.values()
                          if t.get("protein_official") and ALPHA2_RE.search(t["protein_official"][:90])
                          and t["protein_official"].startswith("MV")]
    mam = {k.split("_")[0]: mammal_exon1(loci[k], ALPHA2_REP[k]) for k in ALPHA2_REP}
    report["mammal_alpha2_transcripts"] = alpha2_tx
    report["mammal_exon1alpha2"] = {k: {kk: vv for kk, vv in v.items() if kk not in ("upstream300",)}
                                    for k, v in mam.items()}
    for k, v in mam.items():
        g = loci[f"{k}_SYNE1" if k == "human" else f"{k}_Syne1"].find_tx(
            GIANT_REF[f"{k}_SYNE1" if k == "human" else f"{k}_Syne1"])
        nums = [i + 1 for i, e in enumerate(g["exons"]) if e == v["second_exon"]]
        v["giant_exon_number_of_first_shared_exon"] = nums[0] if nums else None
        print(f"{k:6s} exon1a2 {v['exon']} ({v['length']} bp) 5'UTR {v['utr5']} coding {v['coding_nt']} nt "
              f"phase {v['end_phase']} donor {v['donor']} pep {v['peptide']} | first shared = giant exon "
              f"{v['giant_exon_number_of_first_shared_exon']}")
    mam_peps = {k: v["peptide"] for k, v in mam.items()}
    human_alpha2 = loci["human_SYNE1"].find_tx(ALPHA2_REP["human_SYNE1"])["protein_official"]

    # 2-4. zebrafish
    zf = {}
    for key in ["zebrafish_syne1a", "zebrafish_syne1b", "zebrafish_syne1a_GRCz12tu", "zebrafish_syne1b_GRCz12tu"]:
        loc = loci[key]
        giant = loc.find_tx(GIANT_REF[key])
        fse = first_shared_exon(loc, giant, human_alpha2)
        scan = scan_intron(loc, fse, giant["strand"], mam_peps, mam)
        zf[key] = dict(first_shared_exon=fse, scan={k: v for k, v in scan.items() if k != "intron_seq"},
                       intron_seq=scan["intron_seq"])
        top = scan["candidates"][0] if scan["candidates"] else None
        print(f"\n{key}: first shared exon = {GIANT_REF[key]} exon {fse['exon_number']} {fse['exon']} "
              f"(upstream supplies {fse['upstream_contribution']} nt of split codon {fse['residue_aa']}{fse['residue']}); "
              f"intron {fse['intron']} ({fse['intron_len']} bp); {scan['n_candidates']} ORF candidates")
        if top:
            print(f"   top candidate: {top['genomic']} coding {top['coding_nt']} nt donor {top['donor_seq']} "
                  f"pep {top['peptide']} scores {top['scores']}")
    # decide: syne1a GRCz11 top candidate
    best = zf["zebrafish_syne1a"]["scan"]["candidates"][0]
    loc = loci["zebrafish_syne1a"]
    cds_a, cds_b = best["genomic"]
    coding = loc.sub(cds_a, cds_b, 1)
    assert coding.startswith("ATG") and len(coding) % 3 == 2
    donor = loc.sub(cds_b + 1, cds_b + 8, 1)

    # 5' boundary (approximate TSS) by homology: exon 1alpha2 5'UTRs share an "AAGGGGA" motif
    # with the zebrafish sequence upstream of the ATG. The zebrafish exon start is placed at the
    # same distance upstream of that motif as in the curated human transcript (MANE Plus
    # Clinical NM_001347702.2). Mouse/rat TSSs lie further upstream, so this is a minimal exon.
    up = loc.sub(cds_a - 200, cds_a - 1, 1)
    motif_pos_z = up.rfind("AAGGGGA")
    offsets = {k: v["seq"][:v["utr5"]].rfind("AAGGGGA") for k, v in mam.items()}
    assert motif_pos_z >= 0 and offsets["human"] >= 0
    exon_start = cds_a - 200 + motif_pos_z - offsets["human"]
    exon5 = [exon_start, cds_b]
    report["zebrafish_exon1alpha2"] = dict(
        assembly="GRCz11", chrom="chr20", strand="+", cds_part=[cds_a, cds_b], coding_nt=len(coding),
        predicted_exon=exon5, predicted_exon_note=(
            "3' end (GC splice donor) and CDS start are defined by the ORF; the 5' end (TSS) is approximate, "
            "placed by homology to the AAGGGGA motif that lies "
            + ", ".join(f"{k} {o} nt" for k, o in offsets.items() if o >= 0)
            + " downstream of the exon start (human used as anchor)"),
        coding_seq=coding, peptide=translate(coding[:-2]), donor=donor, donor_type=donor[:2],
        kozak=loc.sub(cds_a - 6, cds_a + 3, 1), exon_seq=loc.sub(exon5[0], exon5[1], 1),
        scores_vs_mammal_peptides=best["scores"])
    print("\nzebrafish syne1a exon1alpha2 (GRCz11): CDS part chr20:%d-%d, predicted exon chr20:%d-%d" %
          (cds_a, cds_b, exon5[0], exon5[1]))

    # 4. whole-gene scans with shuffled controls (fish alpha2 peptide from Betta syne1a, an independent query)
    probe = "MVLDLGSDNMAGCSSEEDLPDCDCDVN"
    report["whole_gene_scans"] = {}
    for key, strand in [("zebrafish_syne1a", 1), ("zebrafish_syne1b", -1),
                        ("zebrafish_syne1a_GRCz12tu", 1), ("zebrafish_syne1b_GRCz12tu", -1)]:
        r = whole_gene_scan(loci[key], strand, probe)
        report["whole_gene_scans"][key] = r
        print(f"whole-gene scan {key}: best {r['best_score']} at {r['best_hit_genomic']} {r['best_hit_seq']} "
              f"(shuffled mean {r['shuffled_mean']}, max {r['shuffled_max']})")
    for k, p in mam_peps.items():
        r = whole_gene_scan(loci["zebrafish_syne1b"], -1, p, n_shuffle=20)
        report["whole_gene_scans"][f"zebrafish_syne1b_vs_{k}"] = r

    # 5. other assemblies (exact match of the coding part + donor)
    probe_nt = coding + donor[:6]
    places = {"GRCz11 (danRer11; Ensembl 116)": dict(acc="chr20", start=cds_a, end=cds_b)}
    l12 = loci["zebrafish_syne1a_GRCz12tu"]
    i = l12.seq.find(probe_nt)
    assert i >= 0
    places["GRCz12tu (NCBI RefSeq)"] = dict(acc="NC_133195.1 (chr20)", start=l12.start + i,
                                            end=l12.start + i + len(coding) - 1)
    off12 = l12.start + i - cds_a
    # GRCz13ab: fetch a window around the syne1a 3' region of NC_141040.1
    w0, w1 = 30200000, 30280876
    s13 = ncbi_region("NC_141040.1", w0, w1)
    j = s13.find(probe_nt)
    if j >= 0:
        places["GRCz13ab (NCBI RefSeq)"] = dict(acc="NC_141040.1 (chr20)", start=w0 + j, end=w0 + j + len(coding) - 1)
    report["zebrafish_exon1alpha2"]["other_assemblies"] = places
    report["zebrafish_exon1alpha2"]["GRCz12tu_predicted_exon"] = [exon5[0] + off12, exon5[1] + off12]
    # check that the predicted exon (incl. 5' part) maps colinearly
    ex12 = l12.sub(exon5[0] + off12, exon5[1] + off12, 1)
    report["zebrafish_exon1alpha2"]["GRCz12tu_exon_identical"] = ex12 == report["zebrafish_exon1alpha2"]["exon_seq"]
    for k, v in places.items():
        print(f"   {k}: {v['acc']}:{v['start']}-{v['end']}")

    # 6. vertebrate survey
    surv = vertebrate_survey()
    report["vertebrate_survey_counts"] = {}
    for r in surv:
        report["vertebrate_survey_counts"][r["group"]] = report["vertebrate_survey_counts"].get(r["group"], 0) + 1
    with open(RESULTS / "vertebrate_alpha2_survey.tsv", "w") as fh:
        fh.write("group\torganism\tgene\taccession\tlength\tunique_N_terminus\tjunction\n")
        for r in surv:
            fh.write(f"{r['group']}\t{r['organism']}\t{r['gene']}\t{r['accession']}\t{r['length']}\t"
                     f"{r['unique_nterm']}\t{r['junction']}\n")
    print("vertebrate survey:", report["vertebrate_survey_counts"])
    report["vertebrate_survey"] = surv

    # zebrafish intron sequences for the viewer
    report["zebrafish_scans"] = {k: dict(first_shared_exon=v["first_shared_exon"], **v["scan"]) for k, v in zf.items()}
    save_json(RESULTS / "alpha2_exon_analysis.json", report)
    write_fasta(RESULTS / "zebrafish_syne1a_intron_hosting_exon1alpha2_GRCz11.fa",
                [("syne1a_intron%d_GRCz11 chr20:%d-%d(+)" % (zf["zebrafish_syne1a"]["first_shared_exon"]["prev_exon_number"],
                                                          *zf["zebrafish_syne1a"]["first_shared_exon"]["intron"]),
                  zf["zebrafish_syne1a"]["intron_seq"])])


if __name__ == "__main__":
    main()
