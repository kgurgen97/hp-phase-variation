#!/usr/bin/env python3
"""G5: run g5caller extract for manifest samples in parallel (evidence caches only; no disease labels are read).
Development samples (CALIBRATION_AND_DEVELOPMENT) are always allowed. POST_FREEZE_PLATFORM_TRANSFER_ONLY samples
are refused unless --frozen-rules PATH is given and its sha256 equals metadata/g5/g5_caller_rules.sha256."""
import argparse, csv, hashlib, os, subprocess, sys
from multiprocessing import Pool
ap = argparse.ArgumentParser()
ap.add_argument('--study', required=True); ap.add_argument('--rules', required=True)
ap.add_argument('--outdir', required=True); ap.add_argument('--frozen-rules'); ap.add_argument('--procs', type=int, default=8)
a = ap.parse_args()
rows = [r for r in csv.DictReader(open('metadata/g5/g5_development_manifest.tsv'), delimiter='\t') if r['bioproject'] == a.study]
assert rows
if rows[0]['g5_use'] == 'POST_FREEZE_PLATFORM_TRANSFER_ONLY':
    assert a.frozen_rules, 'post-freeze samples require --frozen-rules'
    want = open('metadata/g5/g5_caller_rules.sha256').read().split()[0]
    got = hashlib.sha256(open(a.frozen_rules, 'rb').read()).hexdigest()
    assert want == got and a.rules == a.frozen_rules, 'caller rules are not the frozen ones'
os.makedirs(a.outdir, exist_ok=True)
def one(r):
    out = f"{a.outdir}/{r['biosample']}.evidence.tsv.gz"
    if os.path.exists(out + '.meta.json'): return r['biosample'], 'cached'
    fs = r['fastq_paths'].split(';')
    mates = [f for f in fs if f.endswith('_1.fastq.gz') or f.endswith('_2.fastq.gz')]
    unp = [f for f in fs if f not in mates]
    cmd = [sys.executable, '-m', 'g5caller', 'extract', '--catalogue', 'metadata/g5/g5_repeat_catalogue.tsv', '--sample', r['biosample'], '--rules', a.rules, '--out', out]
    if r['library_layout'] == 'SINGLE': cmd += ['--se'] + fs
    else:
        cmd += ['--pe'] + mates
        if unp: cmd += ['--orphan'] + unp
    env = dict(os.environ, PYTHONPATH='scripts/g5')
    p = subprocess.run(cmd, capture_output=True, text=True, env=env)
    return r['biosample'], 'ok' if p.returncode == 0 else 'ERR ' + p.stderr[-300:]
if __name__ == '__main__':
    with Pool(a.procs) as p:
        for b, s in p.imap_unordered(one, rows): print(b, s, flush=True)
