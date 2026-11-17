#!/bin/bash
# M5: rsID-bearing abstracts from the full PubMed 2026 baseline (echo's verified copy), 14 processes
set -euo pipefail
~/Documents/Bioinformatics/bin/python3 "$IN_DIR/m5/01_extract_candidates.py" ~/echo/raw/baseline ~/varclaims/work/cand 14
ls ~/varclaims/work/cand/cand_*.parquet | wc -l > "$OUT_DIR/n_files.txt"
