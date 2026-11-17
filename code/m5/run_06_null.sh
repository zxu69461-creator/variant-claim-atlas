#!/bin/bash
# M5: charter amendment 3 null sample — LD proxies (same settings as run_04_ld.sh) for 20,000 random GWAS Catalog rsIDs (seed 20261107),
# and 1000G EUR and EAS allele frequencies for every claim and null rsID
set -euo pipefail
cd ~/varclaims; mkdir -p work/null; cp "$IN_DIR/null/null_rsids.txt" work/null/
PL=~/Documents/Bioinformatics/bin/plink19
run() { pop=$1; chr=$2; ref=$([ $pop = eur ] && echo ~/bjdt_ref/g1000_eur || echo ~/varclaims/ref/g1000_eas)
  $PL --bfile $ref --chr $chr --r2 --ld-snp-list work/null/null_rsids.txt --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0.8 \
      --threads 1 --memory 1500 --out work/null/${pop}_chr$chr > /dev/null 2>&1 || echo "no result $pop chr$chr"; }
export -f run; export PL
for pop in eur eas; do for c in $(seq 1 22); do echo $pop $c; done; done | xargs -P 14 -n 2 bash -c 'run "$@"' _
for pop in eur eas; do cat work/null/${pop}_chr*.ld | awk 'NR==1 || $1!="CHR_A"' | gzip > "$OUT_DIR/null_${pop}_proxies.ld.gz"; done
cat work/ld/query_rsids.txt work/null/null_rsids.txt | sort -u > work/null/freq_rsids.txt
for pop in eur eas; do ref=$([ $pop = eur ] && echo ~/bjdt_ref/g1000_eur || echo ~/varclaims/ref/g1000_eas)
  $PL --bfile $ref --extract work/null/freq_rsids.txt --freq --threads 4 --out work/null/${pop}_freq > /dev/null
  gzip -c work/null/${pop}_freq.frq > "$OUT_DIR/${pop}_freq.frq.gz"; done
ls -la "$OUT_DIR"
