#!/usr/bin/env python3
"""G5 disease-blind Phase 1 tables from the consolidated locus-level calls (metadata: bioproject/instrument/layout/depth only).
usage: g5_phase1_tables.py LOCUS_CALLS.tsv OUTDIR [BENCH_PER_LOCUS.tsv]"""
import csv, sys, collections, statistics as st, math, os
calls_f, out = sys.argv[1], sys.argv[2]
bench_f = sys.argv[3] if len(sys.argv) > 3 else None
rd = lambda f: list(csv.DictReader(open(f), delimiter='\t'))
calls = rd(calls_f)
cat = rd('metadata/g5/g5_repeat_catalogue.tsv')
qc = {r['biosample']: r for r in rd('metadata/g4/g4_phase1_qc.tsv')}   # technical QC only (no disease columns exist in this file)
A = 25
CATL = collections.defaultdict(list)
for r in cat: CATL[r['locus_id']].append(r)
ELIG = {l for l, v in CATL.items() if any(x['caller_eligible'].startswith('YES') for x in v)}
def rep(l):  # representative member: prefer 26695 then first eligible
    v = [x for x in CATL[l] if x['caller_eligible'].startswith('YES')]
    return sorted(v, key=lambda x: (x['reference_short'] != '26695', x['member_id']))[0]
studies = ['PRJNA622860', 'PRJNA360417', 'PRJNA1103397']
samples = collections.defaultdict(list)
for r in calls: samples[r['bioproject']].append(r)
bs_by_study = {s: sorted({r['biosample'] for r in calls if r['bioproject'] == s}) for s in studies}
def is_call(r): return r['callable'] == 'YES'
# phase1 calls
cols = list(calls[0])
with open(out + '/g5_phase1_calls.tsv', 'w', newline='') as f:
    w = csv.DictWriter(f, cols, delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(calls)
# callability
bysample = collections.defaultdict(list)
for r in calls: bysample[r['biosample']].append(r)
rows = []
for b, rs in sorted(bysample.items()):
    c = [r for r in rs if is_call(r)]
    q = qc[b]
    rows.append(dict(record_type='SAMPLE', bioproject=rs[0]['bioproject'], biosample=b, run_accession=rs[0]['run_accession'], instrument=rs[0]['instrument'], layout=rs[0]['layout'],
                     est_depth_x=q['est_depth_x'], read_len_mean=q['read_len_mean'], eligible_loci=len(rs), callable_loci=len(c), callable_fraction=round(len(c) / len(rs), 4),
                     callable_high=sum(r['confidence_class'] == 'HIGH' for r in c), callable_medium=sum(r['confidence_class'] == 'MEDIUM' for r in c), callable_low=sum(r['confidence_class'] == 'LOW' for r in c),
                     median_spanning_reads=st.median([int(r['spanning_reads']) for r in c]) if c else 0,
                     mixed_flag_MIXED=sum(r['mixed_flag'] == 'MIXED' for r in c), mixture_candidate=sum(r['mixed_flag'] == 'MIXTURE_CANDIDATE' for r in c),
                     known_or_strong_callable=sum(1 for r in c if rep(r['locus_id'])['evidence_class'] in ('KNOWN_PV', 'STRONG_PV_CANDIDATE'))))
for s in studies:
    rr = [r for r in rows if r['bioproject'] == s]
    rows.append(dict(record_type='STUDY_SUMMARY', bioproject=s, biosample=str(len(rr)) + ' samples', instrument=';'.join(sorted({r['instrument'] for r in rr})), layout=';'.join(sorted({r['layout'] for r in rr})),
                     est_depth_x=st.median(float(r['est_depth_x']) for r in rr), read_len_mean=st.median(float(r['read_len_mean']) for r in rr), eligible_loci=len(ELIG),
                     callable_loci=st.median(r['callable_loci'] for r in rr), callable_fraction=round(st.median(r['callable_fraction'] for r in rr), 4),
                     callable_high=st.median(r['callable_high'] for r in rr), callable_medium=st.median(r['callable_medium'] for r in rr), callable_low=st.median(r['callable_low'] for r in rr),
                     median_spanning_reads=st.median(r['median_spanning_reads'] for r in rr), mixed_flag_MIXED=st.median(r['mixed_flag_MIXED'] for r in rr),
                     mixture_candidate=st.median(r['mixture_candidate'] for r in rr), known_or_strong_callable=st.median(r['known_or_strong_callable'] for r in rr)))
c1 = list(rows[0]); 
with open(out + '/g5_phase1_callability.tsv', 'w', newline='') as f:
    w = csv.DictWriter(f, c1, delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(rows)
# platform transfer
def auc(pos, neg):
    if not pos or not neg: return None
    s = 0
    for p in pos:
        for n in neg: s += 1 if p > n else 0.5 if p == n else 0
    return s / (len(pos) * len(neg))
depth = {b: float(qc[b]['est_depth_x']) for b in qc}
rlen = {s: st.median(float(qc[b]['read_len_mean']) for b in bs_by_study[s]) for s in studies}
sdepth = {s: st.median(depth[b] for b in bs_by_study[s]) for s in studies}
MINSP = 8
bloc = collections.defaultdict(list)
for r in calls: bloc[r['locus_id']].append(r)
bench = {}
if bench_f:
    bl = rd(bench_f)
    g = collections.defaultdict(lambda: collections.defaultdict(list))
    for r in bl:
        if r['mix'] == '100:0': g[r['member_id'].rsplit('_', 1)[0]][(r['layout'], int(r['depth']))].append(r['callable'] == 'YES')
    for l, d in g.items(): bench[l] = {k: sum(v) / len(v) for k, v in d.items()}
plat = []
for l in sorted(ELIG):
    m = rep(l); T = int(m['tract_length_bp'])
    rs = bloc.get(l, [])
    d = dict(locus_id=l, evidence_class=m['evidence_class'], context=m['context'], motif_canonical=m['motif_canonical'], motif_length=m['motif_length'], tract_length_bp=T, required_span_bp=T + 2 * A)
    fr = {}
    for s in studies:
        x = [r for r in rs if r['bioproject'] == s]
        c = [r for r in x if is_call(r)]
        fr[s] = len(c) / len(x) if x else None
        d['callable_fraction_' + s] = round(fr[s], 4) if fr[s] is not None else 'NA'
        d['median_spanning_reads_' + s] = st.median([int(r['spanning_reads']) for r in c]) if c else 0
        L = rlen[s]; D = sdepth[s]
        d['expected_spanning_fragments_analytic_' + s] = round(max(0.0, D * (L - T - 2 * A + 1) / L), 1)
    a, b = fr['PRJNA360417'], fr['PRJNA1103397']
    if a is None or b is None: cls = 'NA'
    elif abs(a - b) >= 0.5: cls = 'PLATFORM_CONFOUNDED_CALLABILITY'
    elif a >= 0.8 and b >= 0.8: cls = 'CALLABLE_BOTH'
    elif a < 0.2 and b < 0.2: cls = 'UNCALLABLE_BOTH'
    else: cls = 'PARTIAL'
    d['callability_class'] = cls
    pos = [depth[r['biosample']] for r in rs if is_call(r)]; neg = [depth[r['biosample']] for r in rs if not is_call(r)]
    au = auc(pos, neg) if len(pos) >= 5 and len(neg) >= 5 else None
    d['depth_auc_callable_vs_not'] = round(au, 3) if au is not None else 'NA'
    d['depth_dependent'] = 'YES' if au is not None and au >= 0.8 else 'NO' if au is not None else 'NA'
    e360, e1103 = d['expected_spanning_fragments_analytic_PRJNA360417'], d['expected_spanning_fragments_analytic_PRJNA1103397']
    lim = (e360 < MINSP and e1103 >= 3 * MINSP) or (e1103 < MINSP and e360 >= 3 * MINSP)
    d['predicted_platform_limited'] = 'YES' if lim else 'NO'
    d['span_exceeds_read_SE150'] = 'YES' if T + 2 * A > 150 else 'NO'; d['span_exceeds_read_PE251'] = 'YES' if T + 2 * A > 251 else 'NO'
    bb = bench.get(l)
    d['synthetic_callable_SE150_30x_pure'] = round(bb.get(('SE150', 30), float('nan')), 2) if bb else 'NA'
    d['synthetic_callable_PE250_30x_pure'] = round(bb.get(('PE250', 30), float('nan')), 2) if bb else 'NA'
    plat.append(d)
with open(out + '/g5_platform_transfer.tsv', 'w', newline='') as f:
    w = csv.DictWriter(f, list(plat[0]), delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(plat)
# within host
man = {r['biosample']: r for r in rd('metadata/g3/g3_within_host_validation_manifest.tsv') if r['study'] == 'PRJNA622860'}
wh = []
for l in sorted(ELIG):
    rs = [r for r in bloc.get(l, []) if r['bioproject'] == 'PRJNA622860']
    c = [r for r in rs if is_call(r)]
    m = rep(l)
    alleles = collections.Counter(r['dominant_bp'] for r in c)
    hi = [r for r in c if r['confidence_class'] in ('HIGH', 'MEDIUM') and r['mixed_flag'] == 'NONE']
    al_hi = collections.Counter(r['dominant_bp'] for r in hi)
    site = collections.defaultdict(set)
    for r in hi: site[man[r['biosample']]['region']].add(r['dominant_bp'])
    supp = [k for k, v in al_hi.items() if v >= 2]
    wh.append(dict(locus_id=l, evidence_class=m['evidence_class'], context=m['context'], motif_canonical=m['motif_canonical'], gene_name=m['gene_name'], product=m['product'][:60],
                   isolates_total=len(rs), isolates_callable=len(c), isolates_confident_pure=len(hi),
                   distinct_dominant_alleles=len(alleles), alleles_confident_pure=';'.join(f'{k}bp:{v}' for k, v in sorted(al_hi.items(), key=lambda x: int(x[0]))),
                   variable_confident=('YES' if len(supp) >= 2 and len(hi) >= 10 else 'NO'),
                   isolates_MIXED=sum(r['mixed_flag'] == 'MIXED' for r in c), isolates_MIXTURE_CANDIDATE=sum(r['mixed_flag'] == 'MIXTURE_CANDIDATE' for r in c),
                   functional_states=';'.join(f'{k}:{v}' for k, v in sorted(collections.Counter(r['functional_state'] for r in c).items())),
                   alleles_by_region='|'.join(f'{k}=>{",".join(sorted(v))}' for k, v in sorted(site.items())),
                   note='technical observation only; within-host variation is not interpreted biologically in G5'))
with open(out + '/g5_within_host_repeat_profiles.tsv', 'w', newline='') as f:
    w = csv.DictWriter(f, list(wh[0]), delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(wh)
print('tables written:', len(calls), 'calls;', len(plat), 'platform loci;', sum(1 for r in plat if r['callability_class'] == 'PLATFORM_CONFOUNDED_CALLABILITY'), 'platform-confounded;', sum(1 for r in wh if r['variable_confident'] == 'YES'), 'within-host variable')
