#!/bin/bash
# M5: LD proxies (r2 >= 0.8, 1 Mb) for every rsID named in a candidate abstract, 1000G EUR and EAS (MAGMA g1000 bfiles), per chromosome in parallel
set -euo pipefail
cd ~/varclaims; mkdir -p work/ld
~/Documents/Bioinformatics/bin/python3 -c "
import duckdb
r = duckdb.sql(\"select distinct unnest(rsids) r from read_parquet('work/items_all.parquet')\").fetchall()
open('work/ld/query_rsids.txt', 'w').write(''.join(f'rs{x[0]}\n' for x in r)); print(len(r), 'query rsIDs')"
PL=~/Documents/Bioinformatics/bin/plink19
run() { pop=$1; chr=$2; ref=$([ $pop = eur ] && echo ~/bjdt_ref/g1000_eur || echo ~/varclaims/ref/g1000_eas)
  $PL --bfile $ref --chr $chr --r2 --ld-snp-list work/ld/query_rsids.txt --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0.8 \
      --threads 1 --memory 1500 --out work/ld/${pop}_chr$chr > /dev/null 2>&1 || echo "no result $pop chr$chr"; }
export -f run; export PL
for pop in eur eas; do for c in $(seq 1 22); do echo $pop $c; done; done | xargs -P 14 -n 2 bash -c 'run "$@"' _
for pop in eur eas; do cat work/ld/${pop}_chr*.ld | awk 'NR==1 || $1!="CHR_A"' | gzip > work/ld/${pop}_proxies.ld.gz; done
ls -la work/ld/*_proxies.ld.gz | tee "$OUT_DIR/ld_files.txt"
for pop in eur eas; do echo $pop $(gzip -dc work/ld/${pop}_proxies.ld.gz | wc -l); done | tee -a "$OUT_DIR/ld_files.txt"
