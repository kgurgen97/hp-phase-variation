#!/usr/bin/env python3
"""G7 follow-up: verify exact frozen PRJNA678459 raw-read membership and integrity from existing local
files only (no network call, no download). Expected membership comes only from already-harmonized
metadata/harmonized/run_linkage.tsv and metadata/g7/g7_targeted_raw_read_download_proposal.tsv.
Frozen scope: metadata/g7_followup/g7f_scope_rules.tsv (F03)."""
import csv, gzip, hashlib, os, sys

RAW = 'data/raw/g7_followup/PRJNA678459'
OUT = 'metadata/g7_followup/g7f_frozen_membership_manifest.tsv'

expected = {}
for r in csv.DictReader(open('metadata/harmonized/run_linkage.tsv'), delimiter='\t'):
    if r['bioproject'] == 'PRJNA678459':
        expected[r['run_accession']] = r

est_bytes_by_biosample = {}
for r in csv.DictReader(open('metadata/g3/g3_reserve_manifest.tsv'), delimiter='\t'):
    if r['study'] == 'PRJNA678459':
        est_bytes_by_biosample[r['biosample']] = r['est_bytes']

prop_biosamples = set()
for r in csv.DictReader(open('metadata/g7/g7_targeted_raw_read_download_proposal.tsv'), delimiter='\t'):
    if r['bioproject'] == 'PRJNA678459':
        prop_biosamples = set(r['biosamples'].split(';'))

local_runs = sorted(os.listdir(RAW)) if os.path.isdir(RAW) else []

rows = []
status_all = 'PASS'
for run in sorted(set(expected) | set(local_runs)):
    exp = expected.get(run)
    row = dict(run_accession=run, biosample=exp['biosample'] if exp else 'NOT_IN_FROZEN_EXPECTATION',
               expected_in_run_linkage=('YES' if exp else 'NO'),
               biosample_in_g7_proposal=('YES' if exp and exp['biosample'] in prop_biosamples else 'NO'),
               local_dir_present=('YES' if run in local_runs else 'NO'))
    d = os.path.join(RAW, run)
    files = sorted(os.listdir(d)) if os.path.isdir(d) else []
    row['local_files'] = ';'.join(files) if files else 'NONE'
    row['n_local_files'] = len(files)
    actual_bytes = 0
    gzip_ok = True
    sha_list = []
    for f in files:
        p = os.path.join(d, f)
        actual_bytes += os.path.getsize(p)
        h = hashlib.sha256()
        try:
            with gzip.open(p, 'rb') as fh:
                while True:
                    chunk = fh.read(1024 * 1024)
                    if not chunk:
                        break
                    h.update(chunk)
            sha_list.append(f + ':' + h.hexdigest()[:16])
        except Exception as e:
            gzip_ok = False
            sha_list.append(f + ':GZIP_FAIL:' + type(e).__name__)
    row['actual_bytes'] = actual_bytes
    exp_bytes = est_bytes_by_biosample.get(exp['biosample']) if exp else None
    row['expected_bytes_est'] = exp_bytes if exp_bytes else 'NA'
    row['gzip_integrity'] = 'OK' if gzip_ok else 'FAIL'
    row['sha256_prefix_per_file'] = ';'.join(sha_list)
    if exp_bytes not in (None, '', 'NA'):
        try:
            row['byte_match'] = 'YES' if int(exp_bytes) == actual_bytes else 'NO'
        except ValueError:
            row['byte_match'] = 'NA'
    else:
        row['byte_match'] = 'NA'
    row['md5_check'] = 'NOT_AVAILABLE_NO_NETWORK_ACCESS_THIS_RUN'
    ok = (row['expected_in_run_linkage'] == 'YES' and row['local_dir_present'] == 'YES'
          and row['n_local_files'] == 2 and row['gzip_integrity'] == 'OK' and row['byte_match'] == 'YES')
    row['run_status'] = 'PASS_SIZE_AND_GZIP_ONLY' if ok else 'FAIL_OR_MISSING'
    if row['run_status'] != 'PASS_SIZE_AND_GZIP_ONLY':
        status_all = 'FAIL_OR_MISSING'
    rows.append(row)

cols = list(rows[0].keys())
with open(OUT, 'w', newline='') as fh:
    w = csv.DictWriter(fh, cols, delimiter='\t', lineterminator='\n')
    w.writeheader()
    w.writerows(rows)

n_expected = len(expected)
n_local = len(local_runs)
n_pass = sum(1 for r in rows if r['run_status'] == 'PASS_SIZE_AND_GZIP_ONLY')
print('expected_runs', n_expected, 'local_run_dirs', n_local, 'pass', n_pass, 'of', len(rows))
print('OVERALL_MEMBERSHIP_STATUS', status_all if n_pass == n_expected == n_local == len(rows) else 'FAIL_OR_MISSING')
