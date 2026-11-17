#!/bin/bash
# Amendment 4 chain, run once after the reader check is submitted:
#   collect answers -> 10 reader check (gate) -> [pass] 08 heldout matching -> 11 heldout confirmatory estimates
# Usage: scripts/run_chain.sh DIR_WITH_ANSWER_JSON   (ArtifactData out_dir dump: <dir>/answers/<item id>.json)
set -uo pipefail
cd /Volumes/ResearchData/Claude/varclaims
PY=/Library/Frameworks/Python.framework/Versions/3.14/bin/python3
LOG=work/audit/chain.log; exec > >(tee -a "$LOG") 2>&1
echo "== chain start $(date '+%F %T')"
$PY - "$1" <<'PYEOF' || { echo "CHAIN_FAIL collect"; exit 1; }
import json, sys, pathlib
d = pathlib.Path(sys.argv[1]); out = {}
for f in sorted(d.rglob('*.json')):
    j = json.load(open(f)); doc = j.get('data', j); out[j.get('id', f.stem)] = doc
items = {i['id'] for i in json.load(open('work/audit/items_v1.json'))}
out = {k: v for k, v in out.items() if k in items}
json.dump(out, open('work/audit/responses_v1.json', 'w'), ensure_ascii=False)
print('collected', len(out), 'of', len(items))
PYEOF
$PY scripts/10_reader_check.py; rc=$?
if [ $rc -eq 2 ]; then echo "CHAIN_STOP incomplete answers"; exit 2; fi
if [ $rc -eq 3 ]; then echo "CHAIN_STOP gate failed: no held-out matching (amendment 4.1)"; exit 3; fi
if [ $rc -ne 0 ]; then echo "CHAIN_FAIL reader check rc=$rc"; exit 1; fi
$PY scripts/08_match.py heldout || { echo "CHAIN_FAIL 08 heldout"; exit 1; }
$PY scripts/11_confirm.py heldout || { echo "CHAIN_FAIL 11 heldout"; exit 1; }
echo "CHAIN_DONE $(date '+%F %T')"
# confirmatory figures and manuscript draft (do not change any estimate; a failure here leaves the estimates in place)
ok=1
for f in 22_fig2_support 23_fig3_size 24_fig4_time 25_figS1_sensitivity; do $PY scripts/$f.py heldout || ok=0; done
$PY scripts/30_build_manuscript.py heldout || ok=0
$PY ~/.claude/skills/manuscript/lint.py manuscript/varclaims_heldout.md | grep -E "FAIL" || true
[ $ok -eq 1 ] && echo "DRAFT_DONE" || echo "DRAFT_FAIL (estimates are complete; fix the figure or text step)"
