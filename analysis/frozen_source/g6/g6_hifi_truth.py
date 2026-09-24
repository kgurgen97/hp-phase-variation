#!/usr/bin/env python3
"""G6 HiFi truth per frozen rules H06 to H10. usage: g6_hifi_truth.py EVIDENCE_DIR OUT.tsv"""
import sys, csv, gzip, glob, collections, os
ev_dir, out = sys.argv[1], sys.argv[2]
cat = {r['member_id']: r for r in csv.DictReader(open('metadata/g5/g5_repeat_catalogue.tsv'), delimiter='\t')}
locmem = collections.defaultdict(list)
for m, r in cat.items(): locmem[r['locus_id']].append(m)
def rep(l):
    v = [cat[m] for m in locmem[l] if cat[m]['caller_eligible'].startswith('YES')]
    return sorted(v, key=lambda x: (x['reference_short'] != '26695', x['member_id']))[0] if v else None
def strand_ok(rd, allele):   # each strand >= 2 reads for that allele
    p = sum(1 for l, s in rd if l == allele and s == '+'); m = sum(1 for l, s in rd if l == allele and s == '-'); return p >= 2 and m >= 2
def truth(rd, k):
    lens = [l for l, _ in rd]; n = len(lens); c = collections.Counter(lens)
    if n < 10: return 'UNCALLABLE_HIFI', [], f'{n} spanning reads (<10)'
    if n < 20: return 'AMBIGUOUS_HIFI', [], f'{n} spanning reads (10 to 19)'
    top, ct = c.most_common(1)[0]
    for a, v in c.most_common()[1:]:
        apart = abs(a - top) >= 2 if k == 1 else (a % k == 0 and top % k == 0 and abs(a - top) >= k)
        if apart and ct / n >= 0.15 and v / n >= 0.15 and ct >= 5 and v >= 5 and strand_ok(rd, top) and strand_ok(rd, a) and (ct + v) / n >= 0.9:
            return 'MIXED_HIFI_TRUTH', sorted([top, a]), f'two alleles {top}:{ct} and {a}:{v} of {n}'
    if ct / n >= 0.9 and (k == 1 or top % k == 0) and strand_ok(rd, top): return 'HIGH_CONFIDENCE_HIFI_TRUTH', [top], f'dominant {top} fraction {ct/n:.3f}'
    why = []
    if ct / n < 0.9: why.append(f'dominant fraction {ct/n:.2f} < 0.90')
    if k > 1 and top % k != 0: why.append('dominant allele not unit-aligned')
    if not strand_ok(rd, top): why.append('dominant allele lacks two reads on each strand')
    return 'AMBIGUOUS_HIFI', [], '; '.join(why) or 'criteria not met'
rows = []
for f in sorted(glob.glob(ev_dir + '/*.tsv.gz')):
    sample = os.path.basename(f).split('.')[0]
    by = collections.defaultdict(lambda: collections.defaultdict(list))
    for l in gzip.open(f, 'rt'):
        e = l.rstrip('\n').split('\t')
        if e[0] == 'sample' or e[3] == '1': continue
        by[e[1]][e[2]].append((int(e[6]), e[5]))
    for loc in sorted(locmem):
        r = rep(loc)
        if r is None: continue
        mem = by.get(loc, {})
        base = dict(isolate=sample, locus_id=loc, motif_length=r['motif_length'], motif_canonical=r['motif_canonical'], tract_length_bp=r['tract_length_bp'], context=r['context'], known_pv='YES' if r['evidence_class'] == 'KNOWN_PV' else 'NO', evidence_class=r['evidence_class'])
        if not mem:
            rows.append(dict(base, member_used='NA', hifi_spanning_depth=0, allele_distribution='NA', dominant_allele_bp='NA', dominant_fraction='NA', strand_plus=0, strand_minus=0, repeat_count='NA', truth_status='UNCALLABLE_HIFI', truth_alleles='NA', truth_confidence='NA', truth_basis='no spanning HiFi reads')); continue
        best = max(mem, key=lambda m: (len(mem[m]), m == loc + '_26695')); rd = mem[best]; k = int(cat[best]['motif_length']); lens = [l for l, _ in rd]
        c = collections.Counter(lens); top, ct = c.most_common(1)[0]
        st, al, why = truth(rd, k)
        for m, v in mem.items():
            if m != best and len(v) >= 0.5 * len(rd) and collections.Counter(x[0] for x in v).most_common(1)[0][0] != top and st in ('HIGH_CONFIDENCE_HIFI_TRUTH', 'MIXED_HIFI_TRUTH'): st, al, why = 'AMBIGUOUS_HIFI', [], 'member allele conflict (preserved)'
        sc = collections.Counter(s for _, s in rd)
        rows.append(dict(base, member_used=best, hifi_spanning_depth=len(rd), allele_distribution=';'.join(f'{a}bp:{v}' for a, v in sorted(c.items())), dominant_allele_bp=top, dominant_fraction=round(ct / len(rd), 3), strand_plus=sc.get('+', 0), strand_minus=sc.get('-', 0),
                         repeat_count=round(top / k, 2), truth_status=st, truth_alleles=','.join(map(str, al)) or 'NA', truth_confidence=('HIGH' if st == 'HIGH_CONFIDENCE_HIFI_TRUTH' else 'MIXED' if st == 'MIXED_HIFI_TRUTH' else 'LOW' if st == 'AMBIGUOUS_HIFI' else 'NA'), truth_basis=why))
with open(out, 'w', newline='') as fh:
    w = csv.DictWriter(fh, list(rows[0]), delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(rows)
print(len(rows), collections.Counter(r['truth_status'] for r in rows))
