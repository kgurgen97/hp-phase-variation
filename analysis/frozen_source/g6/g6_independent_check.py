#!/usr/bin/env python3
"""G6 independent boundary check (rules I01 to I09, metadata/g6/g6_independent_check_rules.tsv). NEW code: does not import g5caller or any earlier
extractor. bwa mem -x pacbio + samtools; allele length from the net indel content of a 200 bp window around the tract (no anchors).
usage: g6_independent_check.py ISOLATE [MAXREADS]   -> results/g6/independent_check/<isolate>.tsv"""
import sys, os, csv, gzip, re, subprocess, collections
iso = sys.argv[1]; maxreads = int(sys.argv[2]) if len(sys.argv) > 2 else 0
out_dir = 'results/g6/independent_check'; os.makedirs(out_dir, exist_ok=True)
rd = lambda f: list(csv.DictReader(open(f), delimiter='\t'))
cat = {r['member_id']: r for r in rd('metadata/g5/g5_repeat_catalogue.tsv')}
chk = [r for r in rd('results/g6/g6_hifi_nonreference_validation.tsv') if r['isolate'] == iso and (r['allele_class'] == 'NONREFERENCE' or r['known_pv'] == 'YES')]
# genome (panel) reader: simple, independent
def fasta(path):
    seqs = {}; n = None
    for l in gzip.open(path, 'rt'):
        if l.startswith('>'): n = l[1:].split()[0]; seqs[n] = []
        else: seqs[n].append(l.strip().upper())
    return {k: ''.join(v) for k, v in seqs.items()}
G = fasta('data/reference/g5_panel/panel_all.fna.gz')
FL, CTX = 100, 250
ref_fa = f'{out_dir}/{iso}.contexts.fa'; info = {}
with open(ref_fa, 'w') as f:
    for r in chk:
        m = cat[r['member_used']]; s = int(m['start']); e = int(m['end']); c0 = max(1, s - CTX); c1 = min(len(G[m['contig']]), e + CTX)
        name = f"{r['locus_id']}_{m['reference_short']}"; info[name] = dict(row=r, t0=s - c0, t1=e - c0 + 1, ref_tract=e - s + 1, ctx_len=c1 - c0 + 1)
        f.write(f'>{name}\n{G[m["contig"]][c0 - 1:c1]}\n')
man = [x for x in rd('metadata/g6/g6_hifi_manifest.tsv') if x['isolate'] == iso and x['read_type'] == 'LONG_HIFI'][0]
subprocess.run(['bwa', 'index', ref_fa], capture_output=True, check=True)
reads = man['local_path']
if maxreads:
    sub = f'{out_dir}/{iso}.subset.fq.gz'; n = 0
    with gzip.open(reads, 'rt') as f, gzip.open(sub, 'wt') as o:
        for i, l in enumerate(f):
            o.write(l)
            if i % 4 == 3:
                n += 1
                if n >= maxreads: break
    reads = sub
p1 = subprocess.Popen(['bwa', 'mem', '-x', 'pacbio', '-t', '2', ref_fa, reads], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
p2 = subprocess.Popen(['samtools', 'view', '-F', '2308', '-q', '20'], stdin=p1.stdout, stdout=subprocess.PIPE, text=True)
CIG = re.compile(r'(\d+)([MIDNSHP=X])')
samfh = gzip.open(f'{out_dir}/{iso}.context_alignments.sam.gz', 'wt')
cov = collections.Counter(); clean = collections.Counter(); dirty = collections.Counter(); alle = collections.defaultdict(collections.Counter)
for line in p2.stdout:
    c = line.rstrip('\n').split('\t'); name = c[2]
    if name not in info: continue
    samfh.write(line)
    I = info[name]; wL = I['t0'] - FL; wR = I['t1'] + FL
    if wL < 0 or wR > I['ctx_len']: continue
    pos = int(c[3]) - 1; ops = [(int(n), o) for n, o in CIG.findall(c[5])]
    rpos = 0; ref = pos; rL = rR = None; fl_indels = 0
    for n, o in ops:
        if o in 'M=X':
            if rL is None and ref <= wL < ref + n: rL = rpos + (wL - ref)
            if rR is None and ref < wR <= ref + n: rR = rpos + (wR - ref)
            ref += n; rpos += n
        elif o == 'I':
            # flank insertion: strictly inside a flank zone; an insertion AT the tract boundary is a tract-region event (aligners left-align indels next to a homopolymer)
            if (wL < ref < I['t0']) or (I['t1'] < ref < wR): fl_indels += 1
            rpos += n
        elif o == 'D':
            # flank deletion: lies entirely inside a flank zone and does not touch the tract
            if (ref >= wL and ref + n < I['t0']) or (ref > I['t1'] and ref + n <= wR): fl_indels += 1
            if rL is None and ref <= wL < ref + n: rL = rpos
            if rR is None and ref < wR <= ref + n: rR = rpos
            ref += n
        elif o == 'S': rpos += n
    if rL is None or rR is None: continue
    cov[name] += 1
    if fl_indels: dirty[name] += 1; continue
    clean[name] += 1; alle[name][I['ref_tract'] + (rR - rL) - (wR - wL)] += 1
samfh.close()
rows = []
for name, I in info.items():
    r = I['row']; c = alle.get(name, collections.Counter()); n = sum(c.values()); top, ct = (c.most_common(1)[0] if c else ('NA', 0)); orig = int(r['hifi_truth_bp'])
    frac = ct / n if n else 0
    if n >= 10 and frac >= 0.8: st = 'CONFIRMED' if top == orig else 'DISAGREES'
    else: st = 'AMBIGUOUS_INDEPENDENT_CHECK'
    rows.append(dict(sample=iso, locus_id=r['locus_id'], motif=f"{r['motif_class']}(k={r['motif_length']})", original_hifi_truth=orig, independent_truth=top, independent_support_reads=n, independent_flank_support=f"{clean.get(name, 0)}/{cov.get(name, 0)}",
                     agreement_status=st, notes=f"dominant fraction {frac:.3f}; distribution {';'.join(f'{k}bp:{v}' for k, v in sorted(c.items()))}; reads excluded for flank indels {dirty.get(name, 0)}; short_call_bp {r['short_dominant_bp']}; known_pv {r['known_pv']}; allele_class {r['allele_class']}; reference_tract_bp {r['reference_tract_bp']}"))
with open(f'{out_dir}/{iso}.tsv', 'w', newline='') as f:
    w = csv.DictWriter(f, list(rows[0]), delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(rows)
print(iso, len(rows), collections.Counter(r['agreement_status'] for r in rows), flush=True)
