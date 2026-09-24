#!/usr/bin/env python3
"""G7 follow-up raw-read QC for PRJNA678459. Reuses the frozen G4 QC rule thresholds unchanged
(metadata/g4/g4_qc_rules.tsv, sha256-verified). See metadata/g7_followup/g7f_scope_rules.tsv F04 for the
disclosed, revised split between full-file fastp metrics and bounded-prefix Python checks (runtime-bounded,
synchronous, no background jobs, per this gate's constraints). Q01/Q06 deviations are documented in F03/F04."""
import csv, gzip, json, os, subprocess
from collections import Counter
from multiprocessing import Pool

OUT = 'results/g7_followup/qc'
M = 'metadata/g7_followup/g7f_qc_manifest.tsv'
K = 21
PREFIX_RECORDS = 200000
FASTP_READS_TO_PROCESS = 500000  # read PAIRS; bounded subsample, see g7f_scope_rules.tsv F04
COMP = str.maketrans('ACGT', 'TGCA')


def canon(s):
    r = s.translate(COMP)[::-1]
    return s if s < r else r


def load_refs():
    S = set()
    d = 'data/reference/g4_identity'
    for fn in sorted(os.listdir(d)):
        if not fn.endswith('.fna.gz'):
            continue
        seq = ''.join(l.strip() for l in gzip.open(os.path.join(d, fn), 'rt') if not l.startswith('>')).upper()
        for i in range(len(seq) - K + 1):
            k = seq[i:i + K]
            if 'N' not in k:
                S.add(canon(k))
    return S


REF = None


def parse_prefix(path, max_records, sample_reads=None, want_names=False):
    st = dict(reads=0, bases=0, n=0, lens=Counter(), err='', names=[], sample=[])
    try:
        with gzip.open(path, 'rt', errors='strict') as f:
            while st['reads'] < max_records:
                h = f.readline()
                if not h:
                    break
                s = f.readline().rstrip('\n')
                p = f.readline()
                q = f.readline().rstrip('\n')
                if not q and not s and not p:
                    break
                if not h.startswith('@') or not p.startswith('+') or len(s) != len(q) or not q:
                    st['err'] = f'BAD_RECORD_{st["reads"] + 1}'
                    break
                st['reads'] += 1
                L = len(s)
                st['bases'] += L
                st['lens'][L] += 1
                st['n'] += s.count('N')
                if want_names:
                    st['names'].append(h.split()[0].rstrip('\n').removesuffix('/1').removesuffix('/2'))
                if sample_reads and len(st['sample']) < sample_reads:
                    st['sample'].append(s)
    except (OSError, EOFError, UnicodeDecodeError) as e:
        st['err'] = 'GZIP_OR_DECODE_ERROR:' + type(e).__name__
    return st


def one(r):
    global REF
    if REF is None:
        REF = load_refs()
    run = r['run_accession']
    pair_files = [r[f'local_path_{i}'] for i in (1, 2) if r[f'local_path_{i}'] != 'NA']
    res = dict(bioproject=r['bioproject'], biosample=r['biosample'], run_accession=run, cohort_role=r['cohort_role'],
               ena_layout=r['library_layout'], instrument=r['instrument'], checksum_status=r['checksum_status'],
               actual_bytes=sum(os.path.getsize(p) for p in pair_files if os.path.exists(p)), fails=[], warns=[])
    res['q01_note'] = 'PASS_SIZE_AND_GZIP_ONLY_NO_MD5_NO_NETWORK'
    if r['checksum_status'] != 'PASS_SIZE_AND_GZIP_ONLY':
        res['fails'].append('Q01')
    if len(pair_files) != 2:
        res['fails'].append('Q03')

    # Q02/Q05/Q10: bounded prefix Python check (disclosed, PARTIAL_PREFIX_CHECK; see g7f_scope_rules.tsv F04)
    ps = [parse_prefix(p, PREFIX_RECORDS, sample_reads=20000 if i == 0 else None, want_names=True) for i, p in enumerate(pair_files)]
    errs = [x['err'] for x in ps if x['err']]
    res['format_errors_prefix'] = ';'.join(errs) or 'NONE'
    res['q02_note'] = f'PARTIAL_PREFIX_CHECK_first_{PREFIX_RECORDS}_pairs'
    if errs or any(x['reads'] == 0 for x in ps):
        res['fails'].append('Q02')
    if ps[0]['reads'] != ps[1]['reads'] or ps[0]['names'] != ps[1]['names']:
        res['fails'].append('Q05')
    res['q05_note'] = f'PARTIAL_PREFIX_CHECK_first_{PREFIX_RECORDS}_pairs'
    prefix_bases = sum(x['bases'] for x in ps)
    prefix_n = sum(x['n'] for x in ps)
    res['n_fraction_prefix'] = round(prefix_n / prefix_bases, 6) if prefix_bases else 0
    res['q10_note'] = f'PARTIAL_PREFIX_CHECK_first_{PREFIX_RECORDS}_pairs'
    if prefix_bases and prefix_n / prefix_bases > 0.05:
        res['fails'].append('Q10')
    elif prefix_bases and prefix_n / prefix_bases > 0.01:
        res['warns'].append('Q10')

    # Q06: ena_read_count not available offline for PRJNA678459 this run; reported NA, not compared (disclosed deviation)
    res['q06_note'] = 'NOT_COMPARED_ENA_READ_COUNT_UNAVAILABLE_OFFLINE'

    # Q14 identity k-mer screen (unchanged: 20000-read sample per the frozen rule)
    hit = 0
    tot = 0
    for s in ps[0]['sample']:
        ks = [s[i:i + K] for i in range(0, len(s) - K + 1)]
        ks = [k for k in ks if 'N' not in k]
        if not ks:
            continue
        tot += 1
        if sum(canon(k) in REF for k in ks) / len(ks) >= 0.30:
            hit += 1
    res['identity_hplike_read_fraction'] = round(hit / tot, 4) if tot else 0
    if res['identity_hplike_read_fraction'] < 0.20:
        res['fails'].append('Q14')
    elif res['identity_hplike_read_fraction'] < 0.60:
        res['warns'].append('Q14')

    # Q07-Q09,Q11-Q13: fastp on a bounded read-pair subsample (full-file fastp timed out; see F04)
    jf = f'{OUT}/fastp/{run}.json'
    os.makedirs(f'{OUT}/fastp', exist_ok=True)
    cmd = ['fastp', '-i', pair_files[0], '-I', pair_files[1], '-j', jf, '-h', f'{OUT}/fastp/{run}.html', '-w', '4', '-Q', '-L', '-G',
           '--reads_to_process', str(FASTP_READS_TO_PROCESS)]
    fp = subprocess.run(cmd, capture_output=True, text=True)
    res['fastp_note'] = f'BOUNDED_SUBSAMPLE_reads_to_process_{FASTP_READS_TO_PROCESS}_pairs_not_full_file'
    if fp.returncode == 0 and os.path.exists(jf):
        j = json.load(open(jf))
        b = j['summary']['before_filtering']
        res['total_reads'] = b['total_reads']
        res['total_bases'] = b['total_bases']
        res['est_depth_x_bounded_subsample'] = round(b['total_bases'] / 1.65e6, 1)
        res['read_len_mean'] = round((b.get('read1_mean_length', 0) + b.get('read2_mean_length', 0)) / 2, 1)
        res['gc_fraction'] = round(b['gc_content'], 4)
        res['q20_fraction'] = round(b['q20_rate'], 4)
        res['q30_fraction'] = round(b['q30_rate'], 4)
        ad = j.get('adapter_cutting', {}).get('adapter_trimmed_reads', 0)
        res['adapter_trimmed_read_fraction'] = round(ad / b['total_reads'], 4) if b['total_reads'] else 0
        res['duplication_rate'] = round(j.get('duplication', {}).get('rate', 0), 4)
        tb = res['total_bases']
        if tb < 8.3e6:
            res['fails'].append('Q07')
        elif tb < 33e6:
            res['warns'].append('Q07')
        if res['total_reads'] and res['read_len_mean'] < 75:
            res['warns'].append('Q08')
        if res['q30_fraction'] < 0.50:
            res['fails'].append('Q09')
        elif res['q30_fraction'] < 0.75:
            res['warns'].append('Q09')
        g = res['gc_fraction']
        if g and (g < 0.30 or g > 0.50):
            res['fails'].append('Q11')
        elif g and (g < 0.33 or g > 0.45):
            res['warns'].append('Q11')
        if res['adapter_trimmed_read_fraction'] > 0.10:
            res['warns'].append('Q12')
        if res['duplication_rate'] > 0.50:
            res['warns'].append('Q13')
    else:
        res['fastp_error'] = fp.stderr[-200:]
        res['fails'].append('Q02')

    res['status'] = 'FAIL' if res['fails'] else ('WARNING' if res['warns'] else 'PASS')
    res['rules_failed'] = ','.join(sorted(set(res['fails']))) or 'NA'
    res['rules_warned'] = ','.join(sorted(set(res['warns']))) or 'NA'
    del res['fails'], res['warns']
    return res


if __name__ == '__main__':
    import hashlib
    want = open('metadata/g4/g4_qc_rules.sha256').read().split()[0]
    got = hashlib.sha256(open('metadata/g4/g4_qc_rules.tsv', 'rb').read()).hexdigest()
    assert want == got, 'frozen G4 QC rules changed; refusing to run'
    os.makedirs(OUT, exist_ok=True)
    rows = [r for r in csv.DictReader(open(M), delimiter='\t')]
    with Pool(6) as p:
        res = p.map(one, rows, chunksize=1)
    cols = list(res[0])
    for x in res:
        for c in x:
            if c not in cols:
                cols.append(c)
    with open('results/g7_followup/qc/g7f_qc.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, cols, delimiter='\t', lineterminator='\n', restval='NA')
        w.writeheader()
        w.writerows(res)
    print(Counter(x['status'] for x in res))
