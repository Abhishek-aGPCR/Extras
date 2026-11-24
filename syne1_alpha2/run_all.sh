#!/usr/bin/env bash
# Rebuild everything: data download (cached in data/cache/), analysis, models, alignments, primers,
# protein domains, c.23560G>T verification, the RNA-seq junction test, the browser, and the variant report (Node).
set -euo pipefail
cd "$(dirname "$0")"
for s in 01_fetch_data 02_find_alpha2_exon 03_build_models 04_alignments 05_primers 06_protein_domains 07_verify_c23560; do
  echo "== $s"
  python3 "scripts/$s.py"
done
# Step 10 (RNA-seq test of exon 1alpha2). The reads it needs are already in data/rnaseq/hits/, so only the
# analysis is re-run here. To re-download them from ENA (~35 GB streamed, ~1.5 h at 6 MB/s), run first:
#   python3 scripts/10_rnaseq_junctions.py patterns && python3 scripts/10_rnaseq_junctions.py scan
echo "== 10_rnaseq_junctions analyse"
python3 scripts/10_rnaseq_junctions.py analyse
echo "== 08_build_html"
python3 scripts/08_build_html.py
echo "== 09_variant_report"
node scripts/09_variant_report.js
