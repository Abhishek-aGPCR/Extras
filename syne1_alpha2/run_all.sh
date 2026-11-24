#!/usr/bin/env bash
# Rebuild everything: data download (cached in data/cache/), analysis, models, alignments, primers,
# protein domains, c.23560G>T verification, the browser, and the variant report (Node).
set -euo pipefail
cd "$(dirname "$0")"
for s in 01_fetch_data 02_find_alpha2_exon 03_build_models 04_alignments 05_primers 06_protein_domains 07_verify_c23560 08_build_html; do
  echo "== $s"
  python3 "scripts/$s.py"
done
echo "== 09_variant_report"
node scripts/09_variant_report.js
