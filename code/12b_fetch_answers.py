"""Pull the reading-check answers of reader slot 1 (the first author, amendment 4.1) from the site database into
DIR/answers/<item id>.json for run_chain.sh. Usage: python 12b_fetch_answers.py DIR [--local]. Prints whether the answers are submitted."""
import json, subprocess, sys
from pathlib import Path
out = Path(sys.argv[1]); local = '--local' in sys.argv
site = Path(__file__).resolve().parents[1] / 'audit_site'
r = subprocess.run(['npx', 'wrangler', 'd1', 'execute', 'DB', '--local' if local else '--remote', '--json', '--command',
                    'SELECT slot, data, updated_at, submitted_at FROM progress WHERE slot = 1'], cwd=site, capture_output=True, text=True, check=True)
rows = json.loads(r.stdout)[0]['results']
if not rows: sys.exit('no saved progress for reader 1')
row = rows[0]; d = json.loads(row['data'])
(out / 'answers').mkdir(parents=True, exist_ok=True)
for k, v in d['answers'].items():
    (out / 'answers' / f'{k}.json').write_text(json.dumps({'id': k, 'data': v}, ensure_ascii=False))
(out / 'state.json').write_text(json.dumps({'updated_at': row['updated_at'], 'submitted_at': row['submitted_at'], 'n': len(d['answers'])}))
print('answers', len(d['answers']), '| submitted', row['submitted_at'])
