"""Candidate abstracts for variant-disease association claims, from the PubMed 2026 baseline (M5 job).
Usage: 01_extract_candidates.py BASELINE_DIR OUTDIR NPROC. Recall-first: a record is kept if title+abstract name at least one dbSNP rsID.
Claim reading (disease, allele, direction, population, significance) happens later. Writes OUTDIR/cand_<seq>.parquet with the rsIDs and
MeSH descriptors kept for mapping."""
import sys, os, gzip, re, glob, time, multiprocessing as mp
import xml.etree.ElementTree as ET
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pmparse as L
RS = re.compile(r"(?i)\brs[ -]?(\d{3,10})\b")
def work(job):
    seq, path, outdir = job
    out = os.path.join(outdir, f"cand_{seq:04d}.parquet")
    if os.path.exists(out): return seq, 'skip', 0
    import pyarrow as pa, pyarrow.parquet as pq
    C = {k: [] for k in ('pmid', 'version', 'pub_year', 'journal', 'issn', 'country', 'ptypes', 'retracted', 'rsids', 'mesh', 'title', 'abstract')}
    cty = L.load_countries(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'countries_v1.tsv'))
    with gzip.open(path) as fh:
        for ev, el in ET.iterparse(fh, events=('end',)):
            if el.tag != 'PubmedArticle': continue
            mc = el.find('MedlineCitation'); a = mc.find('Article'); pm = mc.find('PMID')
            title = L.txt(a.find('ArticleTitle')); ab = ' '.join(L.txt(t) for t in a.findall('Abstract/AbstractText'))
            rs = sorted({int(m) for m in RS.findall(title + ' ' + ab)})
            if rs and ab:
                mesh = [d.text for d in mc.findall('MeshHeadingList/MeshHeading/DescriptorName') if d.text]
                ptypes = '|'.join(t.text for t in a.findall('PublicationTypeList/PublicationType') if t.text)
                auths = a.findall('AuthorList/Author'); aff = L.txt(auths[0].find('AffiliationInfo/Affiliation')) if auths else ''
                ret = int(any(c.get('RefType') == 'RetractionIn' for c in mc.findall('CommentsCorrectionsList/CommentsCorrections')))
                j = a.find('Journal')
                for k, v in zip(C, (int(pm.text), int(pm.get('Version') or 1), L.pub_year(a), j.findtext('Title') or '', j.findtext('ISSN') or '',
                                    L.country(aff, cty), ptypes, ret, rs, '|'.join(mesh), title, ab)):
                    C[k].append(v)
            el.clear()
    pq.write_table(pa.table(C), out + '.tmp'); os.replace(out + '.tmp', out)
    return seq, 'ok', len(C['pmid'])
if __name__ == '__main__':
    base, outdir, nproc = sys.argv[1], sys.argv[2], int(sys.argv[3]); os.makedirs(outdir, exist_ok=True)
    files = sorted(glob.glob(os.path.join(base, 'pubmed26n*.xml.gz')))
    jobs = [(int(re.search(r'n(\d{4})', f).group(1)), f, outdir) for f in files]
    t0 = time.time(); n = 0
    with mp.Pool(nproc) as p:
        for seq, st, k in p.imap_unordered(work, jobs):
            n += k; print(seq, st, k, f'{time.time() - t0:.0f}s', flush=True)
    print('done candidates', n)
