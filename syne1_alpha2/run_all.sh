#!/usr/bin/env bash
# Rebuild everything: data download (cached in data/cache/), analysis, models, alignments, primers, viewer.
set -euo pipefail
cd "$(dirname "$0")"
for s in 01_fetch_data 02_find_alpha2_exon 03_build_models 04_alignments 05_primers 06_build_html; do
  echo "== $s"
  python3 "scripts/$s.py"
done
