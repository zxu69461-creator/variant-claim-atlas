#!/bin/bash
# M5: rsID candidates + PubTator3 disease tags, single process
set -euo pipefail
~/Documents/Bioinformatics/bin/python3 "$IN_DIR/m5/02_items.py" ~/varclaims/work/cand ~/echo/raw/pubtator3/disease2pubtator3.gz ~/varclaims/work/items_all.parquet | tee "$OUT_DIR/items_log.txt"
