#!/bin/bash
# M5: varclaims reader v2 (compact JSON, salvage parse, 2,500 tokens): tuning sample again, then every candidate (development split first). Resumable.
set -euo pipefail
cd ~/varclaims
R=~/depclaims/venv_mlx/bin/python; M=~/models/Qwen3.5-9B-MLX-4bit
$R "$IN_DIR/m5/05_llm_read.py" $M work/tune200.parquet work/tune200_v2.jsonl 8
cp work/tune200_v2.jsonl "$OUT_DIR/"
~/Documents/Bioinformatics/bin/python3 - <<'PY'
import pandas as pd, hashlib
it = pd.read_parquet('work/items_all.parquet')
it['dev'] = it.pmid.map(lambda p: hashlib.md5(str(p).encode()).hexdigest()[-1] in '01')
it.sort_values(['dev', 'pmid'], ascending=[False, True])[['pmid', 'title', 'abstract', 'rsids', 'disease_names', 'dev']].to_parquet('work/items_reader.parquet', index=False)
print(len(it), int(it.dev.sum()))
PY
$R "$IN_DIR/m5/05_llm_read.py" $M work/items_reader.parquet work/reader_all_v2.jsonl 8
