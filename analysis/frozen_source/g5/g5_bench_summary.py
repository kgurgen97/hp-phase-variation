#!/usr/bin/env python3
"""Summaries of the synthetic benchmark per-locus table. usage: g5_bench_summary.py PER_LOCUS.tsv OUT_PREFIX"""
import csv, sys, collections, json
per = list(csv.DictReader(open(sys.argv[1]), delimiter='\t')); out = sys.argv[2]
def stats(rs):
    n = len(rs); c = [r for r in rs if r['callable'] == 'YES']
    ok = [r for r in c if r['dominant_correct'] == 'YES']
    conf = [r for r in c if r['confidence'] in ('HIGH', 'MEDIUM')]
    fcc = [r for r in conf if r['dominant_correct'] == 'NO']
    st = [r for r in c if r['true_state_major'] in ('ON', 'OFF')]
    stok = [r for r in st if r['called_state'] == r['true_state_major']]
    stwrong = [r for r in st if r['called_state'] in ('ON', 'OFF') and r['called_state'] != r['true_state_major']]
    return dict(n=n, callable=len(c), uncallable_rate=round(1 - len(c) / n, 4) if n else 0,
                exact_allele_accuracy=round(len(ok) / len(c), 4) if c else 'NA',
                false_confident_call_rate=round(len(fcc) / len(conf), 4) if conf else 'NA', n_confident=len(conf),
                onoff_truth_callable=len(st), onoff_state_correct=len(stok), onoff_state_wrong_confident=len(stwrong),
                onoff_state_accuracy=round(len(stok) / len(st), 4) if st else 'NA')
def write(fn, groups):
    rows = []
    for k, rs in sorted(groups.items()):
        d = stats(rs); d.update(dict(zip(hdr, k))); rows.append(d)
    cols = list(hdr) + [c for c in rows[0] if c not in hdr]
    with open(fn, 'w', newline='') as f:
        w = csv.DictWriter(f, cols, delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(rows)
    return rows
pure = [r for r in per if r['mix'] == '100:0']
hdr = ('stratum', 'layout', 'depth'); g = collections.defaultdict(list)
for r in pure: g[(r['stratum'], r['layout'], int(r['depth']))].append(r)
write(out + '_pure_by_stratum_layout_depth.tsv', g)
hdr = ('layout', 'depth', 'mix'); g = collections.defaultdict(list)
for r in per: g[(r['layout'], int(r['depth']), r['mix'])].append(r)
write(out + '_by_layout_depth_mix.tsv', g)
hdr = ('stratum',); g = collections.defaultdict(list)
for r in per: g[(r['stratum'],)].append(r)
write(out + '_by_stratum_all.tsv', g)
# mixed performance
M = collections.defaultdict(lambda: collections.Counter())
for r in per:
    if r['callable'] != 'YES': M[(r['mix'], r['depth'])]['uncallable'] += 1; continue
    M[(r['mix'], r['depth'])]['callable'] += 1
    M[(r['mix'], r['depth'])][r['mixed_flag']] += 1
with open(out + '_mixed_performance.tsv', 'w') as f:
    f.write('mix\tdepth\tcallable\tuncallable\tflag_MIXED\tflag_MIXTURE_CANDIDATE\tflag_none\n')
    for (m, d), c in sorted(M.items(), key=lambda x: (x[0][0], int(x[0][1]))):
        f.write('\t'.join(map(str, [m, d, c['callable'], c['uncallable'], c['MIXED'], c['MIXTURE_CANDIDATE'], sum(v for k, v in c.items() if k not in ('callable', 'uncallable', 'MIXED', 'MIXTURE_CANDIDATE'))])) + '\n')
tot = stats(per); tp = stats(pure)
json.dump(dict(all=tot, pure=tp), open(out + '_overall.json', 'w'), indent=1)
print(json.dumps(dict(all=tot, pure=tp)))
