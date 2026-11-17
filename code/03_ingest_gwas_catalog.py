"""Ingest the GWAS Catalog full association table (release 2026-10-29, ontology-annotated alt-full) as Parquet + DuckDB view, and an
rsID-level table (one row per association x mapped EFO trait) for the arbiter. Also ingests the EFO trait-to-parent mapping."""
import sys, subprocess, io
sys.path.insert(0, '/Volumes/ResearchData/Claude/rdata')
import pandas as pd
from rdata.proj import Project
P = Project(__file__); R = P.root / 'raw'
with P.run(__file__, params={'release': '2026-10-29', 'file': 'gwas-catalog-download-associations-alt-full.tsv'}) as run:
    z = R / 'gwas-catalog-associations_ontology-annotated-full.zip'
    run.input('GWAS_Catalog_202611', path=z, url='https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations_ontology-annotated-full.zip')
    raw = subprocess.run(['unzip', '-p', str(z)], capture_output=True, check=True).stdout
    a = pd.read_csv(io.BytesIO(raw), sep='\t', dtype=str, low_memory=False)
    a.columns = [c.lower().replace(' ', '_').replace('/', '_').replace('(', '').replace(')', '').replace('-', '_').replace('%', 'pct')
                 .replace('[', '').replace(']', '') for c in a.columns]
    run.parquet('gwas_catalog_assoc', a, note='GWAS Catalog full associations, release 2026-10-29, all columns as text')
    s = a[['study_accession', 'pubmedid', 'date', 'disease_trait', 'initial_sample_size', 'snp_id_current', 'snps',
           'strongest_snp_risk_allele', 'pvalue_mlog', 'or_or_beta', 'mapped_trait', 'mapped_trait_uri', 'chr_id', 'chr_pos']].copy()
    s['pvalue_mlog'] = pd.to_numeric(s.pvalue_mlog, errors='coerce')
    s['rsid'] = pd.to_numeric(s.snp_id_current, errors='coerce')
    s = s.assign(efo_uri=s.mapped_trait_uri.str.split(', ')).explode('efo_uri')  # trait names contain commas; names come from gwas_efo_parent
    s = s[s.rsid.notna()].copy(); s['rsid'] = s.rsid.astype('int64')
    run.table('gwas_catalog_rsid_trait', s.drop(columns=['mapped_trait_uri']),
              note='GWAS Catalog associations with a current rsID, one row per association x mapped EFO trait')
    m = pd.read_csv(R / 'gwas-efo-trait-mappings.tsv', sep='\t', dtype=str)
    m.columns = ['disease_trait', 'efo_term', 'efo_uri', 'parent_term', 'parent_uri']
    run.table('gwas_efo_parent', m, note='GWAS Catalog reported trait -> EFO term -> EFO parent category')
    print(len(a), 'associations;', len(s), 'rsID x trait rows;', s.rsid.nunique(), 'rsIDs;', s.efo_uri.nunique(), 'EFO terms')
