#!/usr/bin/env python3
"""G7 follow-up: apply the frozen G5 spanning-read caller (scripts/g5/g5caller, rules unchanged, sha256-verified)
to the 10 verified local PRJNA678459 samples. Extraction bounded to 500000 read fragments (pairs) per sample
for synchronous, non-backgrounded runtime (g7f_scope_rules.tsv F05). Restricts reporting to the 261
PRIMARY_TECHNICAL + 39 PROVISIONAL_TECHNICAL loci from metadata/g7/g7_locus_analysis_universe.tsv.
No disease label is read anywhere in this script (F06)."""
import csv, hashlib, os, sys, time
from multiprocessing import Pool

sys.path.insert(0, 'scripts/g5')
from g5caller.__main__ import main as g5main  # noqa: E402

CATALOGUE = 'metadata/g5/g5_repeat_catalogue.tsv'
RULES = 'metadata/g5/g5_caller_rules.tsv'
REFERENCE = 'data/reference/g5_panel/GCF_000008525.1/GCF_000008525.1_ASM852v1_genomic.fna.gz'
RAW = 'data/raw/g7_followup/PRJNA678459'
WORK = 'results/g7_followup/work'
MAXFRAG = 500000


def check_frozen():
    want = open('metadata/g5/g5_caller_rules.sha256').read().split()[0]
    got = hashlib.sha256(open(RULES, 'rb').read()).hexdigest()
    assert want == got, 'frozen G5 caller rules changed; refusing to run'


def one(row):
    run, biosample = row['run_accession'], row['biosample']
    r1 = os.path.join(RAW, run, run + '_1.fastq.gz')
    r2 = os.path.join(RAW, run, run + '_2.fastq.gz')
    ev = f'{WORK}/evidence/{biosample}.evidence.tsv.gz'
    calls = f'{WORK}/calls/{biosample}.calls.tsv'
    os.makedirs(f'{WORK}/evidence', exist_ok=True)
    os.makedirs(f'{WORK}/calls', exist_ok=True)
    t0 = time.time()
    g5main(['extract', '--catalogue', CATALOGUE, '--sample', biosample, '--pe', r1, r2,
            '--rules', RULES, '--out', ev, '--max-fragments', str(MAXFRAG)])
    g5main(['call', '--catalogue', CATALOGUE, '--evidence', ev, '--reference', REFERENCE,
            '--rules', RULES, '--out', calls])
    return biosample, run, round(time.time() - t0, 1)


if __name__ == '__main__':
    check_frozen()
    rows = [r for r in csv.DictReader(open('metadata/g7_followup/g7f_qc_manifest.tsv'), delimiter='\t')]
    with Pool(5) as p:
        for biosample, run, secs in p.imap_unordered(one, rows):
            print(biosample, run, secs, 's', flush=True)
