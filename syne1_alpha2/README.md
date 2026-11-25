# SYNE1 / nesprin-1α2 in zebrafish

Where is the muscle-specific nesprin-1α2 ("SYNE1α2") isoform in zebrafish, and what does it look like next to
the human, mouse and rat genes? This folder holds the scripts, the data they produce and an interactive browser
(`SYNE1_isoform_browser.html`).

## Result

**Nesprin-1α2 is encoded by *syne1a* (chr20). *syne1b* (chr17) has no α2 exon.**

In mammals, nesprin-1α2 starts from an internal promoter. Its unique first exon (exon 1α2) lies in the intron
before human exon 131 and is spliced to the last 16 exons of the giant transcript. The zebrafish ortholog of
exon 1α2 sits in the homologous *syne1a* intron, between Ensembl syne1a-201 exons 127 and 128. Neither Ensembl
nor NCBI RefSeq annotates it.

| Assembly | Exon 1α2 (predicted) | Coding part (86 nt) | Strand |
|---|---|---|---|
| **GRCz11** (UCSC danRer11, Ensembl 116) | **chr20:26,341,126-26,341,308** | chr20:26,341,223-26,341,308 | + |
| GRCz12tu (current NCBI RefSeq) | NC_133195.1:29,055,664-29,055,846 | 29,055,761-29,055,846 | + |
| GRCz13ab (NCBI RefSeq, AB strain) | – | NC_141040.1:30,257,700-30,257,785 | + |

The start codon and the 3′ end (a GC splice donor, `CAG|GCAGGA`) are fixed by the reading frame. The 5′ end
(transcription start) is approximate; it is placed by homology to human exon 1α2.

Exon 1α2, GRCz11 (lowercase = predicted 5′UTR, uppercase = coding):

```
gcaacaccggacactgttatccagactaaagagaaggggaatttctacagggtggaagtttttctttattctggcacactcgattggatacggctcaATGGTGCTGG
ATATGGGCACTGATAACATGGCAGGATGTCCTGGACAGGACTCAGATCTGCCCGAGTGCGACTGTGACATTAACAG|gcaggaaaa...
```

Unique N-terminus: `MVLDMGTDNMAGCPGQDSDLPECDCDIN|RVKKLRETLVAVQQ…` (human: `MVVAEDLSALRMAEDSCVDADLPDCNCDVT|RVKKLKETLVAVQQ…`)

The predicted zebrafish α2 protein is **956 aa** with the DV23 exon (the form expected in muscle) and 933 aa
without it. It is 62.9% identical and 76.8% similar to human nesprin-1α2 (NP_001334631.1, 975 aa). Full
sequences (mRNA, CDS, protein, ± DV23) are in `results/zebrafish_syne1a_alpha2_predicted.fa`. Transcript models
are in `results/*_GRCz11.gtf`, `.bed` (a UCSC custom track) and the GRCz12tu equivalents.

### Evidence

* **Same intron, same frame.** Exon 1α2 has 92 coding nt in human, mouse and rat and 86 nt in zebrafish. All
  four end in phase 2, completing the split codon of the first shared exon. That exon is 163 bp in both species
  and also starts in phase 2: human exon 131, zebrafish syne1a exon 128.
* **Conserved peptide.** This ORF is the best-scoring of 139 ATG-to-donor ORFs in the intron (Smith-Waterman
  56/64/62 against human/mouse/rat). A scan of the whole 268-kb *syne1a* locus with a fish α2 peptide hits this
  exon with a score of 115. The next-best hit scores 38; shuffled-peptide controls average 38 (maximum 47).
* **Conserved non-canonical GC donor.** Human, mouse, rat and zebrafish all splice exon 1α2 with `AG|GCAGG`.
  Gene predictors expect GT donors, which may explain why the exon was missed.
* **Other fish agree.** RefSeq has α2 models for pike, Betta and sterlet *syne1a* and for herring. Their N-termini
  (`MVLDLGTDKMAGSPSQGSPLPECDCDVN…` in herring) closely match the zebrafish prediction. Across vertebrates, 205
  RefSeq proteins of the α2 type were found: 116 mammals, 56 birds, 22 amphibians, 7 reptiles, 4 ray-finned fish.
* **Identical in three assemblies** (GRCz11, GRCz12tu, GRCz13ab).
* ***syne1b* has no candidate.** Its homologous intron and whole locus score only at background level.
* **Expressed in RNA-seq — in the heart** (see below). This is direct transcript evidence, not just a prediction.

### Transcript evidence (RNA-seq)

The predicted exon is a **used exon**. Streaming public RNA-seq from ENA and counting reads that span the
exon 1α2 → exon 128 junction (minimap2 `splice:sr`, duplicate-removed fragments, confirmed by an alignment-free
exact 24-mer; `scripts/10_rnaseq_junctions.py`) shows nesprin-1α2 is expressed in zebrafish, and specifically in
the **heart**:

| Zebrafish tissue (PRJNA1257549, WT adult) | α2 junction (24-mer) | giant junction | exon 1α2 coverage | α2 share |
|---|---|---|---|---|
| fast skeletal muscle (4 fish) | 0 (0) | 40 | 0× | 0% |
| slow skeletal muscle (3 fish) | 2 (2) | 30 | 0.6× | 6% |
| **ventricle / cardiac (4 fish)** | **72 (67)** | 89 | **19×** | **45%** |

All four ventricle fish are positive (per-fish α2 share 19–58 %); the reads join exon 1α2 to exon 128 exactly at
the predicted GC-AG donor (`…TGCGACTGTGACATTAACAG | GGTGAAGAAGTT…`), and the DV23 muscle exon is fully included
(PSI 1.0). The exon is near-absent from bulk adult skeletal muscle — the tissue pattern of cardiac nesprin-1α2 in
mammals, where every cardiomyocyte expresses it but adult skeletal muscle confines it to neuromuscular-junction
nuclei (Holt et al. 2019). Looking only at adult skeletal muscle would have missed it.

**Cross-species and developmental controls** (the same pipeline made strand-generic, `scripts/11_rnaseq_crossspecies.py`):
the method detects α2 where it is independently known, and the minus-strand human/mouse loci count identically to
the plus-strand zebrafish locus.

| Panel | α2 junction (24-mer) | giant | α2 share | role |
|---|---|---|---|---|
| Human healthy left ventricle (×3) | 19 (23) | 35 | 35 % | positive control — detected |
| Mouse C2C12 **myoblast** | **0 (0)** | 5 | 0 % | undifferentiated — absent |
| Mouse C2C12 **myotube, day 3** | **47 (48)** | 15 | 76 % | differentiated — induced |
| Zebrafish 4 dpf whole larva (×4) | 7 (8) | 91 | 7 % | larval — low but detected, exon covered 7× |

The mouse **myoblast → myotube** step (0 → 47) catches the developmental induction of nesprin-1α2 directly, as
expected (Stroud et al. 2017), and doubles as a positive control. The conserved GC-AG α2 donor is recovered in
human, mouse and zebrafish alike by the splice-aware aligner. In whole 4 dpf wild-type larvae the α2 junction is
low but detected (7 fragments, 24-mer 8), where the adult **skeletal** muscle is essentially null — consistent
with α2 appearing during active larval myogenesis (as in the mouse myotubes) and then becoming cardiac-restricted
in the adult. A larval heart dissection or a staged myogenic time-course would resolve the exact timing.

### akap6 and nesprin-1α2 across zebrafish tissues — a rewiring from mammals

In mammals AKAP6/mAKAP is the nuclear-envelope scaffold that nesprin-1α recruits in **cardiac** myocytes
(Pare et al. 2005), and the akap6–nesprin-1α2 module is found in striated muscle and osteoclasts. We measured
**akap6** (cDNA 25-mer capture, RPKM-like; `scripts/12_akap6_expression.py`) next to the **nesprin-1α2** junction
across zebrafish tissues:

| | fast muscle | slow muscle | ventricle (heart) | osteoclast (Ctsk⁺) |
|---|---|---|---|---|
| **akap6** (RPKM) | 1.1 | 1.7 | **0.1** | 0.6 (noisy) |
| **nesprin-1α2** (α2 junction) | 0 | 2 | **72** | 0 |
| *syne1a total* (giant junction) | 40 | 30 | 89 | **2** (≈absent) |

Tissue identities are confirmed by markers (ventricle myl7 = 14,431, mylpfa = 0; osteoclast ctsk = 463). The two
partners of the mammalian complex sit in **different zebrafish tissues**: **akap6 is skeletal-muscle-biased and
~14× lower in the heart** (0.1 vs 1.1–1.7, at the neuronal-background level) — i.e. **not expressed in the
zebrafish heart** — while **nesprin-1α2 is cardiac-specific** (72 in ventricle, ~0 in skeletal muscle). In
osteoclasts both are low and *syne1a itself is barely expressed*, so α2 there is "not detected, low sensitivity".
So the mammalian cardiac nesprin-1α–mAKAP signalosome is **not co-assembled in any zebrafish tissue examined**:
in the heart nesprin-1α2 is **orphaned** (no akap6 partner), and in skeletal muscle akap6 must partner a different
nesprin-1 isoform. (One library per tissue, first pass; akap6 is a low-abundance scaffold, so the heart≪muscle
*ratio* is the robust readout. Re-run with `--full` for replicates.)

Confirm at the bench by RT-PCR with the intron-spanning, cDNA-specific primers in
`results/zebrafish_alpha2_primer_suggestions.tsv` (check specificity with Primer-BLAST first); 5′RACE from exon
128 would pin the transcription start. For example:

| Forward (exon 1α2) | Reverse | cDNA product |
|---|---|---|
| GCTCAATGGTGCTGGATATG | TTCTCTGGATCTCGTGTGAG (exon 2) | 237 bp |
| CTGTGACATTAACAGGGTGAAG (exon 1α2/exon 2 junction) | TTCTCTGGATCTCGTGTGAG | 161 bp |

## Human SYNE1 c.23560G>T (p.Glu7854*)

This is the homozygous nonsense variant of the congenital muscular dystrophy patient myoblast line used by
Holt et al. 2016 and Gimpel et al. 2017. The full report is `results/variant_c23560G-T/report.md`, with figures
`figure_human_c23560G-T.png` and `figure_zebrafish_syne1a_equivalent.png`.

**Position (verified).**

| Coordinate | Value |
|---|---|
| GRCh38 | NC_000006.12:g.152176461C>A (chr6:152,176,461 C>A) |
| GRCh37 | chr6:152,497,596 C>A |
| Transcript | NM_182961.4:c.23560G>T, p.Glu7854Ter, exon 130 of 146 |
| Position in exon | base 100 of 167, 67 nt before the exon end, so no splice-site involvement |
| Codon | GAG>TAG |

It was checked against NM_182961.3 and .4 (identical CDS), the GRCh38 genome, Ensembl variant recoder and VEP,
the Ensembl assembly mapper and ClinVar.

| Isoform | Contains the variant? | Predicted mRNA fate |
|---|---|---|
| Giant (NM_182961.4 and all giant variants) | yes, p.Glu7854* in spectrin repeat 68 | PTC ~2.6 kb upstream of the last exon junction → **NMD** |
| **Nesprin-1α2** (NM_001347702.2, Ensembl SYNE1-202/236) | **no**: 8,149 bp upstream of the α2 transcription start | mRNA and protein sequence unchanged |
| RefSeq isoform 3 (NM_001347701.2) | yes, p.Glu56* | NMD |
| Variant 2 (NM_033071.5) | yes, p.Glu7783* | NMD |
| All 81 annotated human transcripts | 55 carry the stop | all 55 NMD-predicted |

**What the literature gets right and wrong.**
- The giant mRNA was measured 78% lower, as NMD predicts.
- α2 mRNA was 85% lower in patient myotubes. The variant cannot cause that through the α2 mRNA, which does not
  contain it. The report lists testable explanations: silencing of the locus by the premature stop
  (speculative), feedback from losing the exon-130 isoforms, or an effect of this single cell line.
- Holt 2016 writes the codon as "gaa to taa". The reference is GAG, so the change is GAG>TAG, probably a typo.
- The variant is not in ClinVar, and the original clinical report could not be retrieved.

**Zebrafish equivalent.** Human E7854 aligns to **syne1a Glu7748**, codon GAA in exon 127 at GRCz11
chr20:26,339,444-26,339,446. This exon sits just upstream of the intron that holds the zebrafish exon 1α2, as in
human. The analogous allele is ENSDART00000139350.4:c.23242G>T, which gives the ochre stop TAA. It is predicted to
send the syne1a giant mRNA to NMD and leave α2 intact (1,682 bp upstream of the α2 start). syne1b has Asp at this
position.

## The browser

Open `SYNE1_isoform_browser.html` in any modern browser. It is one self-contained file and works offline;
fonts load from Google Fonts when online.

* **Gene and assembly:** human SYNE1 (GRCh38), mouse Syne1 (GRCm39), rat Syne1 (GRCr8), zebrafish syne1a and
  syne1b on GRCz11 (Ensembl, older RefSeq placement) and on GRCz12tu (current RefSeq).
* **Isoforms:** grouped as nesprin-1α2 (muscle), giant, other short KASH, N-terminal and other. Tick an isoform
  to draw it; click a name to focus it.
* **Gene structure:** zoomable exon-intron map (drag to pan, Ctrl/⌘ + scroll or pinch to zoom). Click an exon
  for its coordinates, phase, encoded residues, sequence and UCSC/Ensembl links.
* **Sequence:** genomic (exons uppercase, introns lowercase, + strand or gene orientation), mRNA, CDS with
  translation, and protein. Includes motif search (IUPAC or regex, both strands), go-to position, and copy or
  save as FASTA.
* **Exons:** exon table of the focused isoform, plus the α2 exons compared across species.
* **Alignments:** α2 proteins (10 vertebrates), α2 codon alignment, exon 1α2 nucleotides, α2 N-termini of 206
  vertebrates, and the giant isoforms.
* **Pairwise aligner:** global, local or semi-global affine-gap alignment of any two sequences in the browser
  (protein, CDS, mRNA, an exon, a genomic range) or pasted ones.
* **Variant effect:** enter a variant for any gene/assembly and isoform as
  - coding `c.` (including intronic `c.123+5G>A` and UTR `c.-12`/`c.*30`)
  - transcript `n.`
  - genomic `g.`, `chr:pos ref>alt` or VCF
  - protein `p.` (one- or three-letter codes)
  - `exon N:position`

  Substitutions, del, dup, ins and delins are supported. The tool:
  - checks the reference base against the genome
  - reports HGVS for every isoform, with exon/intron location, codon and amino-acid change, and splice effects
    (including the exon-skipping outcome)
  - predicts NMD for premature stops
  - lists the protein domains lost
  - **for mutations inside the muscle isoform's own exons**, names the functional element hit (KASH anchor,
    LEWD/kinesin motif, transmembrane, a spectrin repeat, or the α2-specific exon-1α2 N-terminus) and what
    losing it does — e.g. a truncation before KASH breaks the LINC complex, a LEWD change abolishes kinesin/MTOC
    recruitment. Pick the muscle isoform in the second menu to enter coordinates relative to it; built-in
    examples cover a KASH missense, a LEWD missense, an NMD-escaping KASH-less truncation, and an exon-1α2
    frameshift.
  - finds the same residue in the other species (one click analyses it there)
  - draws a gene / mRNA-fate / protein-domain figure you can save as SVG or PNG.

  The same engine (`viewer/variant_engine.js`) writes the command-line report.
* **Zebrafish α2 evidence:** coordinates, figure, evidence, primers, methods, data versions and references.

## Reproduce

```bash
pip install -r requirements.txt
```

```bash
bash run_all.sh
```

The steps run in order, and downloads are cached in `data/cache/`:

| Script | Does |
|---|---|
| `scripts/01_fetch_data.py` | Ensembl, NCBI Datasets/E-utilities and UCSC: gene models, loci (gene ± 2 kb), proteins |
| `scripts/02_find_alpha2_exon.py` | mammalian exon 1α2; zebrafish homologous-intron scan; controls; other assemblies; vertebrate survey |
| `scripts/03_build_models.py` | predicted zebrafish α2 models (± DV23), isoform classes, FASTA/GTF/BED |
| `scripts/04_alignments.py` | FAMSA protein alignments, codon and exon 1α2 alignments, exon correspondence |
| `scripts/05_primers.py` | RT-PCR primer suggestions |
| `scripts/06_protein_domains.py` | UniProt Q8NF91 domains (CH, 74 spectrin repeats, LEWD, TM, KASH) projected onto every species |
| `scripts/07_verify_c23560.py` | independent verification of c.23560G>T (RefSeq versions, genome, Ensembl VEP/recoder, GRCh37, ClinVar) |
| `scripts/10_rnaseq_junctions.py` | zebrafish exon 1α2 → exon 128 junction counts in adult muscle/heart RNA-seq (ENA) |
| `scripts/11_rnaseq_crossspecies.py` | the same method across human heart, mouse C2C12, zebrafish larva and osteoclasts (positive/developmental controls) |
| `scripts/12_akap6_expression.py` | zebrafish akap6 vs nesprin-1α2 across fast/slow muscle, ventricle and osteoclasts (RPKM-like cDNA capture) |
| `scripts/08_build_html.py` | builds `SYNE1_isoform_browser.html` from `viewer/template.html` + `viewer/variant_engine.js` |
| `scripts/09_variant_report.js` | Node: variant report, isoform tables and SVG/PNG figures (same engine as the browser) |

Steps 10–11 stream FASTQ from `ftp.sra.ebi.ac.uk` (ENA) and need `mappy` (minimap2) plus `ripgrep`; `run_all.sh`
re-runs only their analysis from the captured reads under `data/rnaseq/` (re-download with the `patterns`/`scan`
subcommands, ~35 GB streamed). Other network access is needed to `rest.ensembl.org`, `eutils.ncbi.nlm.nih.gov`,
`api.ncbi.nlm.nih.gov`, `api.genome.ucsc.edu` and `rest.uniprot.org`. Step 9 needs Node.js; PNG export uses Playwright with Chromium if
installed. Data used here: Ensembl release 116, NCBI RefSeq (downloaded 2026-11-20 and 2026-11-24).

## Files

* `results/alpha2_exon_coordinates.tsv`: exon 1α2 in human, mouse, rat and zebrafish (coordinates, UTR/coding
  lengths, donor, peptide).
* `results/alpha2_exon_analysis.json`: full analysis output (candidate ORFs, scores, controls, survey).
* `results/zebrafish_syne1a_alpha2_predicted.fa`, `*.gtf`, `*.bed`: predicted zebrafish α2 transcripts.
* `results/zebrafish_syne1a_intron_hosting_exon1alpha2_GRCz11.fa`: the scanned intron.
* `results/vertebrate_alpha2_survey.tsv`: α2-type RefSeq proteins across vertebrates.
* `results/alignments/*.fa`: alignments shown in the browser.
* `results/variant_c23560G-T/`: the c.23560G>T files:
  - `report.md` and `report.json`
  - `verification.json`
  - `isoform_consequences_human.tsv` and `isoform_consequences_zebrafish_syne1a.tsv`
  - figures (SVG/PNG)
* `data/domains.json`: protein domain map per species.
* `data/models/`, `data/genomic/`: downloaded gene models and genomic loci; `data/models_built.json`: classified
  models with proteins.

## References

* Holt et al. 2016, BMC Cell Biol 17:26: nesprin-1α2 has a unique 31-aa N-terminus encoded by a sequence
  preceding exon 131. https://doi.org/10.1186/s12860-016-0105-9
* Stroud et al. 2017, J Cell Biol 216:1915: deleting exon 1α2 in mouse causes perinatal lethality and myonuclear
  mispositioning. https://doi.org/10.1083/jcb.201612128
* Duong et al. 2014, PLoS One 9:e94380: α2 is restricted to skeletal and cardiac muscle; DV23 is included in
  muscle. https://doi.org/10.1371/journal.pone.0094380
* Rajgor et al. 2012, PLoS One 7:e40098: nesprin-1 exon numbering and variants. https://doi.org/10.1371/journal.pone.0040098
* Holt et al. 2019, Sci Rep 9:14202. https://doi.org/10.1038/s41598-019-50728-6
* Zhou et al. 2017, Hum Mol Genet 26:2258. https://doi.org/10.1093/hmg/ddx116
* Gimpel et al. 2017, Curr Biol 27:2999: the c.23560G>T patient myotubes lose nesprin-1 at the nuclear envelope
  and microtubule nucleation from the nucleus. https://doi.org/10.1016/j.cub.2017.08.031
* Attali et al. 2009, Hum Mol Genet 18:3462 (a different SYNE1 arthrogryposis mutation). https://doi.org/10.1093/hmg/ddp290
* Baumann et al. 2017, Eur J Hum Genet 25:262: a last-exon SYNE1 nonsense variant keeps wild-type mRNA levels,
  consistent with NMD escape. https://doi.org/10.1038/ejhg.2016.144
* Lindeboom, Supek & Lehner 2016, Nat Genet 48:1112: rules of nonsense-mediated decay. https://doi.org/10.1038/ng.3664
* Lindeboom et al. 2019, Nat Genet 51:1645 (NMDetective): the rule thresholds used here. https://doi.org/10.1038/s41588-019-0517-5
* Bühler et al. 2005, Mol Cell 18:307: PTC-induced transcriptional silencing of Ig minigenes. https://doi.org/10.1016/j.molcel.2005.03.030
* El-Brolosy et al. 2019, Nature 568:193: transcriptional adaptation after mutant mRNA decay (zebrafish, mouse). https://doi.org/10.1038/s41586-019-1064-z
