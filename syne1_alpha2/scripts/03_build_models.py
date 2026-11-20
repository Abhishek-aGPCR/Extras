#!/usr/bin/env python3
"""Step 3 - build the predicted zebrafish syne1a alpha2 transcripts, classify every isoform,
derive mRNA/CDS/protein sequences from the genome and write result files.

Predicted zebrafish alpha2 = exon 1alpha2 (step 2) + the exons of the syne1a giant model from
the first shared exon (homologue of human exon 131) to the last exon. Two variants are made:
with the DV23 exon (the form expected in muscle; DV23 is ~100% included in human muscle,
Duong et al. 2014) and without it.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DATA, RESULTS, load_json, revcomp, save_json, translate, write_fasta  # noqa: E402

MODELS = DATA / "models"
GENOMIC = DATA / "genomic"
ALPHA2_RE = re.compile(r"C[DN]CD[A-Z]{1,2}R[VI]KKLKE")
GIANT_REF = {
    "human_SYNE1": "NM_182961.4", "mouse_Syne1": "NM_001347732.2", "rat_Syne1": "NM_001419849.1",
    "zebrafish_syne1a": "ENSDART00000139350", "zebrafish_syne1b": "ENSDART00000148520",
    "zebrafish_syne1a_GRCz12tu": "XM_073933275.1", "zebrafish_syne1b_GRCz12tu": "XM_087227815.1",
}
DV23_ZF = "TDVVIPESPEAYIKLTENTLRSS"   # zebrafish syne1a DV23 peptide (homologue of human SDVMIPESPEAYVKLTENAIKNT)


def load(key):
    m = load_json(MODELS / f"{key}.json")
    seq = "".join(open(GENOMIC / f"{key}.fa").read().split("\n")[1:])
    return m, seq


def sub(m, seq, a, b, strand=1):
    s = seq[a - m["region"]["start"]:b - m["region"]["start"] + 1]
    return s if strand == 1 else revcomp(s)


def derive(m, seq, t):
    st = t["strand"]
    mrna = "".join(sub(m, seq, s, e, st) for s, e in t["exons"])
    cds = ""
    if t.get("cds"):
        a, b = t["cds"]
        cds = "".join(sub(m, seq, max(s, a), min(e, b), st) for s, e in t["exons"] if max(s, a) <= min(e, b))
    return mrna, cds


def find_tx(m, tid):
    for t in m["transcripts"]:
        if t["id"].split(".")[0] == tid.split(".")[0]:
            return t
    raise KeyError(tid)


def classify(t, giant_prot):
    p = t.get("protein") or ""
    if t.get("predicted"):
        return "alpha2"
    if p.startswith("M") and ALPHA2_RE.search(p[:90]):
        return "alpha2"
    if not p:
        return "other"
    if len(p) >= 7000:
        return "giant"
    tail_same = giant_prot and p[-40:] == giant_prot[-40:]
    head_same = giant_prot and p[:40] == giant_prot[:40]
    if tail_same:
        return "short_kash"
    if head_same:
        return "nterm"
    return "other"


CATEGORY_LABEL = {
    "alpha2": "Nesprin-1α2 (muscle-specific short isoform)",
    "giant": "Nesprin-1 giant (full length)",
    "short_kash": "Other C-terminal (KASH+) short isoforms",
    "nterm": "N-terminal (KASH-less) isoforms",
    "other": "Other / internal isoforms",
}


def protein_exon_map(t):
    """For each exon: [first_aa, last_aa] (1-based) of the protein it encodes, else None."""
    if not t.get("cds"):
        return [None] * len(t["exons"])
    a, b = t["cds"]
    out, acc = [], 0
    for s, e in t["exons"]:
        cs, ce = max(s, a), min(e, b)
        if cs <= ce:
            n = ce - cs + 1
            out.append([acc // 3 + 1, (acc + n - 1) // 3 + 1, acc % 3])
            acc += n
        else:
            out.append(None)
    return out


def main() -> None:
    ana = load_json(RESULTS / "alpha2_exon_analysis.json")
    zx = ana["zebrafish_exon1alpha2"]
    scans = ana["zebrafish_scans"]
    out_views = {}
    fasta_out = []
    gtf_lines = {"GRCz11": [], "GRCz12tu": []}
    bed_lines = {"GRCz11": [], "GRCz12tu": []}

    for key in GIANT_REF:
        m, seq = load(key)
        giant = find_tx(m, GIANT_REF[key])
        # ---- predicted alpha2 models for zebrafish syne1a --------------------------------
        if key in ("zebrafish_syne1a", "zebrafish_syne1a_GRCz12tu"):
            asm = "GRCz11" if key == "zebrafish_syne1a" else "GRCz12tu"
            fse = scans[key]["first_shared_exon"]
            k0 = fse["exon_number"] - 1                      # index of first shared exon
            ex1 = zx["predicted_exon"] if asm == "GRCz11" else zx["GRCz12tu_predicted_exon"]
            cds_start = ex1[1] - zx["coding_nt"] + 1
            downstream = giant["exons"][k0:]
            # locate DV23 exon in the giant model (exon whose CDS encodes the DV23 peptide)
            gp = giant["protein_official"]
            pos = gp.find(DV23_ZF)
            pem = protein_exon_map(giant)
            dv_idx = None
            if pos >= 0:
                for i, r in enumerate(pem):
                    if r and r[0] <= pos + 2 <= r[1] and (giant["exons"][i][1] - giant["exons"][i][0] + 1) == 69:
                        dv_idx = i
            variants = [("plusDV23", downstream)]
            if dv_idx is not None and dv_idx >= k0:
                variants.append(("minusDV23", [e for j, e in enumerate(giant["exons"][k0:], k0) if j != dv_idx]))
            for tag, dex in variants:
                t = dict(id=f"syne1a-alpha2-pred-{tag}-{asm}", source="Predicted (this study)",
                         name=f"syne1a α2 predicted ({'+' if tag == 'plusDV23' else '−'}DV23)",
                         strand=1, exons=[ex1] + dex, cds=[cds_start, giant["cds"][1]],
                         protein_id=None, protein_len=None, predicted=True,
                         tags=["predicted", "muscle-specific α2 model"], biotype="protein_coding",
                         note=("Exon 1α2 identified by homology (this study) + syne1a exons "
                               f"{fse['exon_number']}-{len(giant['exons'])} of {giant['id']}. "
                               "5' end of exon 1α2 approximate."))
                mrna, cds = derive(m, seq, t)
                prot = translate(cds)
                assert prot.endswith("*") and "*" not in prot[:-1], (tag, prot[-20:])
                t["protein"] = prot[:-1]
                t["protein_len"] = len(prot) - 1
                m["transcripts"].insert(0, t)
                if asm == "GRCz11" or True:
                    fasta_out += [(f"{t['id']}_mRNA {asm} chr20 (+) exons={len(t['exons'])}", mrna),
                                  (f"{t['id']}_CDS {asm}", cds),
                                  (f"{t['id']}_protein length={len(prot) - 1}", prot[:-1])]
                chrom = "chr20"
                for i, (s, e) in enumerate(t["exons"], 1):
                    gtf_lines[asm].append(f'{chrom}\tthis_study\texon\t{s}\t{e}\t.\t+\t.\tgene_id "syne1a"; '
                                          f'transcript_id "{t["id"]}"; exon_number "{i}";')
                    cs, ce = max(s, t["cds"][0]), min(e, t["cds"][1])
                    if cs <= ce:
                        gtf_lines[asm].append(f'{chrom}\tthis_study\tCDS\t{cs}\t{ce}\t.\t+\t.\tgene_id "syne1a"; '
                                              f'transcript_id "{t["id"]}"; exon_number "{i}";')
                bs = ",".join(str(e - s + 1) for s, e in t["exons"])
                bo = ",".join(str(s - t["exons"][0][0]) for s, e in t["exons"])
                bed_lines[asm].append(f"{chrom}\t{t['exons'][0][0] - 1}\t{t['exons'][-1][1]}\t{t['id']}\t0\t+\t"
                                      f"{t['cds'][0] - 1}\t{t['cds'][1]}\t200,0,90\t{len(t['exons'])}\t{bs}\t{bo}")
        # ---- every transcript: derived sequences, protein, category, exon->protein map ----
        giant_prot = giant.get("protein_official")
        for t in m["transcripts"]:
            mrna, cds = derive(m, seq, t)
            tp = translate(cds).rstrip("*") if cds else ""
            off = t.get("protein_official")
            if not t.get("predicted"):
                t["protein"] = off or tp
            t["genome_translation_matches_official"] = (off == tp) if off else None
            if off and off != tp:
                t["n_aa_diff_vs_genome"] = sum(1 for a, b in zip(off, tp) if a != b) + abs(len(off) - len(tp))
            t["category"] = classify(t, giant_prot)
            t["exon_aa"] = protein_exon_map(t)
            t["mrna_len"] = len(mrna)
            t.pop("protein_official", None)
        m["giant_ref"] = giant["id"]
        out_views[key] = m
        cats = {}
        for t in m["transcripts"]:
            cats[t["category"]] = cats.get(t["category"], 0) + 1
        mism = [t["id"] for t in m["transcripts"] if t.get("genome_translation_matches_official") is False]
        print(f"{key}: {len(m['transcripts'])} transcripts {cats}; protein!=genome translation: {mism[:6]}")

    # ---- result files ----------------------------------------------------------------
    write_fasta(RESULTS / "zebrafish_syne1a_alpha2_predicted.fa", fasta_out)
    for asm in gtf_lines:
        (RESULTS / f"zebrafish_syne1a_alpha2_predicted_{asm}.gtf").write_text("\n".join(gtf_lines[asm]) + "\n")
        hdr = (f'track name="syne1a_alpha2_pred" description="zebrafish syne1a nesprin-1 alpha2 '
               f'(predicted, {asm})" itemRgb=On\n')
        (RESULTS / f"zebrafish_syne1a_alpha2_predicted_{asm}.bed").write_text(hdr + "\n".join(bed_lines[asm]) + "\n")
    save_json(DATA / "models_built.json", {k: dict(view=v["view"], region=v["region"], giant_ref=v["giant_ref"],
                                                   transcripts=v["transcripts"]) for k, v in out_views.items()})

    # mammalian alpha2 + exon1alpha2 table
    rows = ["species\tassembly\tgene\ttranscript\tprotein\texon1alpha2\tstrand\texon_len\tUTR5_in_exon\t"
            "coding_nt\tsplice_donor\tunique_N_terminus"]
    for sp, v in ana["mammal_exon1alpha2"].items():
        asm = {"human": "GRCh38/hg38 chr6", "mouse": "GRCm39/mm39 chr10", "rat": "GRCr8/rn8 chr1"}[sp]
        rows.append(f"{sp}\t{asm}\t{'SYNE1' if sp == 'human' else 'Syne1'}\t{v['transcript']}\t{v['protein_id']}\t"
                    f"{v['exon'][0]}-{v['exon'][1]}\t-\t{v['length']}\t{v['utr5']}\t{v['coding_nt']}\t"
                    f"{v['donor'][:2]} ({v['donor']})\t{v['peptide']}")
    for asm, (a, b) in [("GRCz11/danRer11 chr20", zx["predicted_exon"]),
                        ("GRCz12tu NC_133195.1", zx["GRCz12tu_predicted_exon"])]:
        rows.append(f"zebrafish\t{asm}\tsyne1a\tpredicted (this study)\t-\t{a}-{b} (CDS part "
                    f"{b - zx['coding_nt'] + 1}-{b})\t+\t{b - a + 1}\t{b - a + 1 - zx['coding_nt']} (approx.)\t"
                    f"{zx['coding_nt']}\t{zx['donor_type']} ({zx['donor']})\t{zx['peptide']}")
    (RESULTS / "alpha2_exon_coordinates.tsv").write_text("\n".join(rows) + "\n")
    print((RESULTS / "alpha2_exon_coordinates.tsv").read_text())


if __name__ == "__main__":
    main()
