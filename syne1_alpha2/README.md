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

### Limitation and how to test it

No zebrafish EST or mRNA (UCSC danRer11) covers the exon, and neither does an Ensembl RNA-seq model. Those models
come only from blastula to 3 dpf embryos and adult male. In mammals α2 is muscle-restricted and switched on
during myotube differentiation, so test it on muscle RNA. `results/zebrafish_alpha2_primer_suggestions.tsv`
lists RT-PCR primer pairs, all intron-spanning and cDNA-specific. Check their specificity with Primer-BLAST
before ordering. For example:

| Forward (exon 1α2) | Reverse | cDNA product |
|---|---|---|
| GCTCAATGGTGCTGGATATG | TTCTCTGGATCTCGTGTGAG (exon 2) | 237 bp |
| CTGTGACATTAACAGGGTGAAG (exon 1α2/exon 2 junction) | TTCTCTGGATCTCGTGTGAG | 161 bp |

5′RACE from exon 128 would define the transcription start.

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
| `scripts/06_build_html.py` | builds `SYNE1_isoform_browser.html` from `viewer/template.html` |

Network access is needed to `rest.ensembl.org`, `eutils.ncbi.nlm.nih.gov`, `api.ncbi.nlm.nih.gov` and
`api.genome.ucsc.edu`. Data used here: Ensembl release 116, NCBI RefSeq (downloaded 2026-11-20).

## Files

* `results/alpha2_exon_coordinates.tsv`: exon 1α2 in human, mouse, rat and zebrafish (coordinates, UTR/coding
  lengths, donor, peptide).
* `results/alpha2_exon_analysis.json`: full analysis output (candidate ORFs, scores, controls, survey).
* `results/zebrafish_syne1a_alpha2_predicted.fa`, `*.gtf`, `*.bed`: predicted zebrafish α2 transcripts.
* `results/zebrafish_syne1a_intron_hosting_exon1alpha2_GRCz11.fa`: the scanned intron.
* `results/vertebrate_alpha2_survey.tsv`: α2-type RefSeq proteins across vertebrates.
* `results/alignments/*.fa`: alignments shown in the browser.
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
