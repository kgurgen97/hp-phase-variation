#!/usr/bin/env python3
"""G5 catalogue v2 (supervisor review): (1) KNOWN_PV only for loci in the reconciliation assignments (no name-only mapping);
(2) STRONG_PV_CANDIDATE only with an explicit comparative-genomic basis (tract-length polymorphism among flank-orthologous panel members);
family-name-only candidates are downgraded to REPEAT_ONLY. Genome-wide discovery thresholds are unchanged. Disease-blind."""
import csv, json, collections
cat = list(csv.DictReader(open('metadata/g5/superseded/g5_repeat_catalogue_v1.tsv'), delimiter='\t'))
assign = json.load(open('metadata/g5/g5_known_pv_assignments.json'))
loci = collections.defaultdict(list)
for r in cat: loci[r['locus_id']].append(r)
ev = collections.Counter(); down = collections.Counter()
for lid, v in loci.items():
    cn = sorted({round(float(x['copy_number']), 2) for x in v})
    anch = any(x['uniquely_anchorable'].startswith('YES') for x in v)
    for r in v:
        old = r['evidence_class']
        ctx_ok = r['context'] in ('CDS', 'PROMOTER_OR_INTERGENIC_UPSTREAM')
        r['known_pv_mapping_basis'] = 'NA'; r['strong_basis'] = 'NA'; r['discovery_status'] = 'GENOME_WIDE_DISCOVERY_RD01_RD10'
        if lid in assign:
            r['evidence_class'] = 'KNOWN_PV'; r['known_pv_mapping_basis'] = assign[lid]['basis']
            r['evidence_source'] = 'literature (reconciled): ' + assign[lid]['entry']
        elif old in ('STRONG_PV_CANDIDATE', 'KNOWN_PV') and len(cn) >= 2 and ctx_ok and r['uniquely_anchorable'].startswith('YES') and r['interruption_flag'] in ('NONE', 'NA') and not r['caller_eligible'].startswith('NO'):
            r['evidence_class'] = 'STRONG_PV_CANDIDATE'
            r['strong_basis'] = f'TRACT_LENGTH_POLYMORPHIC_ACROSS_{len(v)}_PANEL_ORTHOLOGS(copies {cn[0]}-{cn[-1]}; {len(cn)} distinct); slipped-strand-compatible pure tract; context {r["context"]}; disease-independent comparative evidence, not experimental PV'
            r['evidence_source'] = 'RD11a (evidence classes v2)'
        else:
            r['evidence_class'] = 'REPEAT_ONLY'
            r['evidence_source'] = ('downgraded from STRONG_PV_CANDIDATE (v1) because the only basis was gene-family name text or no cross-genome tract-length variation (evidence classes v2)' if old == 'STRONG_PV_CANDIDATE' else
                                    ('KNOWN_PV (v1) mapping not established without a literature-stated tag (evidence classes v2); ' if old == 'KNOWN_PV' else '') + r['evidence_source'])
        if old != r['evidence_class']: down[(old, r['evidence_class'])] += 1
cols = list(cat[0])
for c in ('known_pv_mapping_basis', 'strong_basis', 'discovery_status'):
    if c not in cols: cols.append(c)
with open('metadata/g5/g5_repeat_catalogue.tsv', 'w', newline='') as f:
    w = csv.DictWriter(f, cols, delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(cat)
def lc(v):
    o = {x['evidence_class'] for x in v}; return 'KNOWN_PV' if 'KNOWN_PV' in o else 'STRONG_PV_CANDIDATE' if 'STRONG_PV_CANDIDATE' in o else 'REPEAT_ONLY'
lcount = collections.Counter(lc(v) for v in loci.values()); mcount = collections.Counter(r['evidence_class'] for r in cat)
el = collections.Counter()
for l, v in loci.items(): el[lc(v)] += any(x['caller_eligible'].startswith('YES') for x in v)
print('loci', dict(lcount), 'members', dict(mcount), 'eligible loci by class', dict(el), 'changes', dict(down))
json.dump(dict(loci=dict(lcount), members=dict(mcount), eligible_loci_by_class=dict(el), class_changes={f'{a}->{b}': n for (a, b), n in down.items()}), open('metadata/g5/g5_merge_summary.json', 'w'), indent=1)
# known loci table
kn = list(csv.DictReader(open('metadata/g5/superseded/g5_known_pv_loci_v1.tsv'), delimiter='\t'))
rec = collections.defaultdict(list)
for r in csv.DictReader(open('metadata/g5/g5_known_pv_reconciliation.tsv'), delimiter='\t'):
    rec[r['entry']].append(r)
for r in kn:
    rs = rec[r['gene_name']]
    r['catalogue_locus_ids_matched'] = ';'.join(sorted({x['canonical_locus_id'] for x in rs if x['canonical_locus_id'] != 'NA'})) or 'NONE'
    r['catalogue_match_basis'] = 'see g5_known_pv_reconciliation.tsv'
    r['reconciliation_outcome'] = ';'.join(sorted({x['outcome_class'] for x in rs}))
with open('metadata/g5/g5_known_pv_loci.tsv', 'w', newline='') as f:
    w = csv.DictWriter(f, list(kn[0]), delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(kn)
