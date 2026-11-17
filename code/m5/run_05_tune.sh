#!/bin/bash
# M5: varclaims reader v1 on 200 development-split abstracts (tuning sample, seed 20261104)
set -euo pipefail
cd ~/varclaims
~/Documents/Bioinformatics/bin/python3 - <<'PY'
import pandas as pd, hashlib
it = pd.read_parquet('work/items_all.parquet')
dev = it[it.pmid.map(lambda p: hashlib.md5(str(p).encode()).hexdigest()[-1] in '01')]
dev[dev.rsids.map(len) <= 15].sample(200, random_state=20261104)[['pmid', 'pub_year', 'title', 'abstract', 'rsids', 'disease_names']].to_parquet('work/tune200.parquet', index=False)
print(len(it), len(dev))
PY
~/depclaims/venv_mlx/bin/python "$IN_DIR/m5/05_llm_read.py" ~/models/Qwen3.5-9B-MLX-4bit work/tune200.parquet work/tune200_v1.jsonl 8
cp work/tune200.parquet work/tune200_v1.jsonl "$OUT_DIR/"
