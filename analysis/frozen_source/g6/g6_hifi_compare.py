#!/usr/bin/env python3
"""G6 HiFi comparison and motif validation status per frozen rules H11 to H13 (metadata/g6/g6_high_accuracy_truth_rules.tsv).
usage: g6_hifi_compare.py TRUTH.tsv OUTDIR"""
import sys, csv, glob, collections, statistics as st, math, json
truth_f, out = sys.argv[1], sys.argv[2]
rd = lambda f: list(csv.DictReader(open(f), delimiter='\t'))
truth = {(r['isolate'], r['locus_id']): r for r in rd(truth_f)}
short = {}
for f in sorted(glob.glob('results/g6/hifi_short_calls/*.calls.tsv')):
    by = collections.defaultdict(list)
    for r in rd(f):
        if r['uncallable_reason'].startswith('NOT_CALLER_ELIGIBLE'): continue
        by[r['locus_id']].append(r)
    for loc, rs in by.items():
        cal = [r for r in rs if r['callable'] == 'YES']; pool = cal or rs
        b = max(pool, key=lambda r: (int(r['spanning_reads'] or 0), r['member_id'] == loc + '_26695')); r = dict(b)
        if cal and any(x['member_id'] != b['member_id'] and x['dominant_bp'] != b['dominant_bp'] and int(x['spanning_reads']) >= 0.5 * int(b['spanning_reads']) for x in cal):
            r['callable'] = 'NO'; r['uncallable_reason'] = 'MEMBER_ALLELE_CONFLICT'
        short[(r['sample_id'], loc)] = r
def mclass(t):
    k = int(t['motif_length'])
    return 'MONO_AT' if k == 1 and t['motif_canonical'] == 'A' else 'MONO_GC' if k == 1 else 'DINUCLEOTIDE' if k == 2 else 'TRINUCLEOTIDE' if k == 3 else 'MOTIF_GE4'
def lbin(bp): return 'LE10' if bp <= 10 else '11_15' if bp <= 15 else 'GE16'
def wilson(k, n, z=1.96):
    if n == 0: return (0.0, 1.0)
    p = k / n; d = 1 + z * z / n; c = p + z * z / (2 * n); a = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)); return ((c - a) / d, (c + a) / d)
rows = []
for (iso, loc), t in sorted(truth.items()):
    s = short.get((iso, loc))
    if s is None: continue
    d = dict(isolate=iso, locus_id=loc, motif_class=mclass(t), known_pv=t['known_pv'], context=t['context'], hifi_truth_status=t['truth_status'], hifi_depth=t['hifi_spanning_depth'], hifi_truth_alleles=t['truth_alleles'], hifi_distribution=t['allele_distribution'],
             short_callable=s['callable'], short_dominant_bp=s['dominant_bp'], short_secondary_bp=s['secondary_bp'], short_spanning_reads=s['spanning_reads'], short_confidence=s['confidence_class'], short_mixed_flag=s['mixed_flag'], comparison_universe='NO', outcome='NOT_COMPARED', delta_bp='NA', delta_repeat_units='NA', tract_bin='NA', within_one_unit='NA')
    if s['callable'] == 'YES' and t['truth_status'] in ('HIGH_CONFIDENCE_HIFI_TRUTH', 'MIXED_HIFI_TRUTH'):
        d['comparison_universe'] = 'YES'; al = [int(x) for x in t['truth_alleles'].split(',')]; sd = int(s['dominant_bp']); k = int(t['motif_length']); d['tract_bin'] = lbin(al[0])
        if t['truth_status'] == 'HIGH_CONFIDENCE_HIFI_TRUTH':
            d['delta_bp'] = sd - al[0]; d['delta_repeat_units'] = round((sd - al[0]) / k, 2); d['outcome'] = 'EXACT_CONCORDANT' if sd == al[0] else 'DISCORDANT'; d['within_one_unit'] = 'YES' if abs(sd - al[0]) <= k else 'NO'
        else: d['outcome'] = 'MIXTURE_AGREE' if (s['mixed_flag'] in ('MIXED', 'MIXTURE_CANDIDATE') and sd in al) else 'MIXTURE_NOT_DETECTED'
    elif t['truth_status'] == 'AMBIGUOUS_HIFI': d['outcome'] = 'NOT_SCORED_HIFI_AMBIGUOUS'
    elif t['truth_status'] == 'UNCALLABLE_HIFI': d['outcome'] = 'NOT_SCORED_HIFI_UNCALLABLE'
    elif s['callable'] != 'YES': d['outcome'] = 'NOT_SCORED_SHORT_UNCALLABLE'
    rows.append(d)
with open(out + '/g6_hifi_short_vs_hifi_comparison.tsv', 'w', newline='') as f:
    w = csv.DictWriter(f, list(rows[0]), delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(rows)
hi = [r for r in rows if r['comparison_universe'] == 'YES' and r['hifi_truth_status'] == 'HIGH_CONFIDENCE_HIFI_TRUTH']
catm = {r['member_id']: r for r in rd('metadata/g5/g5_repeat_catalogue.tsv')}
truthm = {(r['isolate'], r['locus_id']): r['member_used'] for r in rd(truth_f)}
def isnonref(r): return int(r['hifi_truth_alleles']) != int(catm[truthm[(r['isolate'], r['locus_id'])]]['tract_length_bp'])
def summ(x):
    ex = [r for r in x if r['outcome'] == 'EXACT_CONCORDANT']; dis = [r for r in x if r['outcome'] == 'DISCORDANT']; lo, hi_ = wilson(len(ex), len(x))
    dl = [int(r['delta_bp']) for r in x]; neg = sum(1 for v in dl if v < 0); pos = sum(1 for v in dl if v > 0)
    return dict(comparisons=len(x), isolates=len({r['isolate'] for r in x}), loci=len({r['locus_id'] for r in x}), exact=len(ex), exact_concordance=round(len(ex) / len(x), 4) if x else 'NA', wilson_low=round(lo, 4), wilson_high=round(hi_, 4),
                within_one_unit=sum(1 for r in x if r['within_one_unit'] == 'YES'), discordant=len(dis), discordant_short_longer=pos, discordant_short_shorter=neg,
                majority_sign_share_of_discordant=round(max(neg, pos) / len(dis), 3) if dis else 'NA', median_delta_bp=st.median(dl) if dl else 'NA', nonreference_allele_comparisons=sum(1 for r in x if isnonref(r)), nonreference_allele_exact=sum(1 for r in x if isnonref(r) and r['outcome'] == 'EXACT_CONCORDANT'), reference_allele_fraction=round(1 - sum(1 for r in x if isnonref(r)) / len(x), 3) if x else 'NA', delta_distribution=';'.join(f'{k}:{v}' for k, v in sorted(collections.Counter(dl).items())))
strata = collections.OrderedDict((k, [r for r in hi if (r['motif_class'] == k)]) for k in ('MONO_AT', 'MONO_GC', 'DINUCLEOTIDE', 'TRINUCLEOTIDE', 'MOTIF_GE4'))
strata['KNOWN_PV'] = [r for r in hi if r['known_pv'] == 'YES']
summary = []
for k, x in strata.items():
    d = summ(x); d['stratum'] = k; summary.append(d)
for b in ('LE10', '11_15', 'GE16'):
    x = [r for r in hi if r['tract_bin'] == b]; d = summ(x); d['stratum'] = 'TRACT_' + b; summary.append(d)
for iso in sorted({r['isolate'] for r in hi}):
    x = [r for r in hi if r['isolate'] == iso]; d = summ(x); d['stratum'] = 'ISOLATE_' + iso; summary.append(d)
d = summ(hi); d['stratum'] = 'ALL_HIGH_CONFIDENCE'; summary.append(d)
cols = ['stratum'] + [c for c in summary[0] if c != 'stratum']
with open(out + '/g6_hifi_summary_by_stratum.tsv', 'w', newline='') as f:
    w = csv.DictWriter(f, cols, delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(summary)
# motif status per H13
NMIN = {'MONO_AT': 100, 'MONO_GC': 50, 'DINUCLEOTIDE': 30, 'TRINUCLEOTIDE': 30, 'MOTIF_GE4': 20, 'KNOWN_PV': 12}
pre = {r['motif_class']: r for r in rd('metadata/g6/g6_motif_validation_status_pre_rescue.tsv')}   # preserved pre-rescue statuses
def status(k, s):
    n = s['comparisons']; nm = NMIN[k]
    if n < 0.5 * nm or (k == 'KNOWN_PV' and s['loci'] < 3 and n < nm): return 'UNVALIDATED_INSUFFICIENT_TRUTH', f'{n} comparisons < 0.5 x {nm}' if n < 0.5 * nm else f'only {s["loci"]} distinct KNOWN_PV loci'
    ex = s['exact_concordance']; lo = s['wilson_low']; dis = s['discordant']; bias = dis >= 20 and s['majority_sign_share_of_discordant'] != 'NA' and s['majority_sign_share_of_discordant'] > 0.7
    if ex < 0.90 or bias: return 'UNVALIDATED_SYSTEMATIC_DISCORDANCE', f'exact {ex:.3f}; discordant {dis}; majority-sign share {s["majority_sign_share_of_discordant"]}'
    if n >= nm and s['isolates'] >= 3 and ex >= 0.95 and lo >= 0.90: return 'VALIDATED_EXACT', f'{n} comparisons, exact {ex:.3f}, Wilson lower {lo:.3f}'
    if n >= nm and s['isolates'] >= 3 and (ex >= 0.90 or lo >= 0.85): return 'VALIDATED_WITH_LIMITATION', f'{n} comparisons, exact {ex:.3f}, Wilson lower {lo:.3f}'
    if n < nm and ex >= 0.95: return 'VALIDATED_WITH_LIMITATION', f'thin evidence: {n} of {nm} required comparisons, exact {ex:.3f}'
    return 'UNVALIDATED_INSUFFICIENT_TRUTH', f'{n} comparisons, exact {ex:.3f}, criteria for validation not met'
res = []
sm = {d['stratum']: d for d in summary}
for k in ('MONO_AT', 'MONO_GC', 'DINUCLEOTIDE', 'TRINUCLEOTIDE', 'MOTIF_GE4', 'KNOWN_PV'):
    s = sm[k]; stt, why = status(k, s); p = pre.get(k, {})
    res.append(dict(motif_class=k, status_before_rescue=p.get('status_before_rescue', 'NA'), basis_before_rescue=p.get('basis', 'NA'), hifi_comparisons=s['comparisons'], hifi_isolates=s['isolates'], hifi_loci=s['loci'], hifi_exact_concordance=s['exact_concordance'], hifi_wilson_low=s['wilson_low'], hifi_discordant=s['discordant'],
                    hifi_median_delta_bp=s['median_delta_bp'], nonreference_allele_comparisons=s['nonreference_allele_comparisons'], nonreference_allele_exact=s['nonreference_allele_exact'], reference_allele_fraction=s['reference_allele_fraction'], hifi_delta_distribution=s['delta_distribution'], n_min_required=NMIN[k], status_after_rescue=stt, status_final_basis=why, interpretation_caveat=(f"{round(100*s['reference_allele_fraction'])}% of comparisons are reference-panel-length alleles (a predict-the-reference baseline would score {s['reference_allele_fraction']}); only {s['nonreference_allele_comparisons']} comparisons tested a non-reference allele ({s['nonreference_allele_exact']} exact). Comparisons cover only loci with resolvable anchors and unambiguous HiFi dominant alleles; both platforms share the catalogue anchor definition; 5 isolates from 5 patients, one region" if s['comparisons'] else 'no comparisons')))
with open('metadata/g6/g6_motif_validation_status.tsv', 'w', newline='') as f:
    w = csv.DictWriter(f, list(res[0]), delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(res)
mx = [r for r in rows if r['comparison_universe'] == 'YES' and r['hifi_truth_status'] == 'MIXED_HIFI_TRUTH']
json.dump(dict(truth_status_counts=dict(collections.Counter(r['hifi_truth_status'] for r in rows)), comparisons_high=len(hi), mixed_truth_comparisons=len(mx), mixture_agree=sum(1 for r in mx if r['outcome'] == 'MIXTURE_AGREE'), motif_status={r['motif_class']: r['status_after_rescue'] for r in res}), open(out + '/g6_hifi_overall.json', 'w'), indent=1)
for r in res: print(r['motif_class'], r['hifi_comparisons'], r['hifi_exact_concordance'], r['status_after_rescue'])
