#!/usr/bin/env python3
"""G7 follow-up: consolidate member-level G5 caller calls to locus-level calls for the 10 PRJNA678459
samples, restricted to the frozen PRIMARY_TECHNICAL (261) and PROVISIONAL_TECHNICAL (39) loci from
metadata/g7/g7_locus_analysis_universe.tsv, reported in SEPARATE tables (kept separate per issue instruction).
No disease label is read anywhere in this script (g7f_scope_rules.tsv F06)."""
import csv, glob, collections, statistics as st

CAT = {r['member_id']: r for r in csv.DictReader(open('metadata/g5/g5_repeat_catalogue.tsv'), delimiter='\t')}
UNIV = {r['locus_id']: r['g7_locus_class'] for r in csv.DictReader(open('metadata/g7/g7_locus_analysis_universe.tsv'), delimiter='\t')}
QC = {r['biosample']: r for r in csv.DictReader(open('results/g7_followup/qc/g7f_qc.tsv'), delimiter='\t')}

locus_members = collections.defaultdict(list)
for m, r in CAT.items():
    locus_members[r['locus_id']].append(m)
locus_eligible = {l: any(CAT[m]['caller_eligible'].startswith('YES') for m in ms) for l, ms in locus_members.items()}

TARGET_TIERS = {'PRIMARY_TECHNICAL', 'PROVISIONAL_TECHNICAL'}
FIELDS = ['biosample', 'run_accession', 'locus_id', 'g7_locus_class', 'member_used', 'member_conflict', 'callable',
          'uncallable_reason', 'spanning_reads', 'dominant_bp', 'dominant_fraction', 'confidence_class', 'mixed_flag',
          'functional_state']

rows_out = []
for f in sorted(glob.glob('results/g7_followup/work/calls/*.calls.tsv')):
    biosample = f.split('/')[-1].replace('.calls.tsv', '')
    run = QC[biosample]['run_accession']
    by = collections.defaultdict(list)
    for r in csv.DictReader(open(f), delimiter='\t'):
        by[r['locus_id']].append(r)
    for loc, rs in by.items():
        tier = UNIV.get(loc)
        if tier not in TARGET_TIERS:
            continue
        if not locus_eligible.get(loc):
            continue
        cal = [r for r in rs if r['callable'] == 'YES']
        pool = cal or rs
        best = max(pool, key=lambda r: (int(r['spanning_reads'] or 0), r['member_id'] == loc + '_26695'))
        conflict = 'NONE'
        if cal:
            oth = [r for r in cal if r['member_id'] != best['member_id'] and r['dominant_bp'] != best['dominant_bp']
                   and int(r['spanning_reads']) >= 0.5 * int(best['spanning_reads'])]
            if oth:
                conflict = 'MEMBER_ALLELE_CONFLICT:' + ','.join(sorted({r['member_id'] for r in oth}))
        o = dict(biosample=biosample, run_accession=run, locus_id=loc, g7_locus_class=tier,
                  member_used=best['member_id'], member_conflict=conflict)
        for k in FIELDS[6:]:
            o[k] = best.get(k, 'NA')
        if conflict != 'NONE':
            o['callable'] = 'NO'
            o['uncallable_reason'] = 'MEMBER_ALLELE_CONFLICT'
            o['functional_state'] = 'UNCALLABLE'
            o['confidence_class'] = 'NA'
        rows_out.append(o)

with open('results/g7_followup/g7f_locus_level_calls.tsv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, FIELDS, delimiter='\t', lineterminator='\n')
    w.writeheader()
    w.writerows(rows_out)

# per-sample callability summary, PRIMARY_TECHNICAL and PROVISIONAL_TECHNICAL kept in separate columns
by_sample = collections.defaultdict(list)
for r in rows_out:
    by_sample[r['biosample']].append(r)
summary = []
for b, rs in sorted(by_sample.items()):
    prim = [r for r in rs if r['g7_locus_class'] == 'PRIMARY_TECHNICAL']
    prov = [r for r in rs if r['g7_locus_class'] == 'PROVISIONAL_TECHNICAL']
    prim_c = [r for r in prim if r['callable'] == 'YES']
    prov_c = [r for r in prov if r['callable'] == 'YES']
    summary.append(dict(
        biosample=b, run_accession=rs[0]['run_accession'],
        bounded_subsample_fragments=500000,
        est_depth_x_bounded_subsample=QC[b].get('est_depth_x_bounded_subsample', 'NA'),
        primary_technical_eligible_loci=len(prim), primary_technical_callable=len(prim_c),
        primary_technical_callable_fraction=round(len(prim_c) / len(prim), 4) if prim else 'NA',
        provisional_technical_eligible_loci=len(prov), provisional_technical_callable=len(prov_c),
        provisional_technical_callable_fraction=round(len(prov_c) / len(prov), 4) if prov else 'NA',
        primary_median_spanning_reads=(st.median(int(r['spanning_reads']) for r in prim_c) if prim_c else 0),
        provisional_median_spanning_reads=(st.median(int(r['spanning_reads']) for r in prov_c) if prov_c else 0),
    ))
cols = list(summary[0].keys())
with open('results/g7_followup/g7f_callability_summary.tsv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, cols, delimiter='\t', lineterminator='\n')
    w.writeheader()
    w.writerows(summary)

print(len(rows_out), 'locus-level rows;', len(summary), 'samples')
print('PRIMARY_TECHNICAL loci in universe:', sum(1 for v in UNIV.values() if v == 'PRIMARY_TECHNICAL'))
print('PROVISIONAL_TECHNICAL loci in universe:', sum(1 for v in UNIV.values() if v == 'PROVISIONAL_TECHNICAL'))
