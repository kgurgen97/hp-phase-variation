#!/usr/bin/env python3
"""G5 locus eligibility (policy metadata/g5/g5_platform_policy.tsv). Disease-blind: uses only the platform audit and caller eligibility."""
import csv, collections
rd = lambda f: list(csv.DictReader(open(f), delimiter='\t'))
cat = rd('metadata/g5/g5_repeat_catalogue.tsv'); plat = {r['locus_id']: r for r in rd('results/g5/g5_platform_transfer.tsv')}
loci = collections.defaultdict(list)
for r in cat: loci[r['locus_id']].append(r)
def cls(v):
    o = {x['evidence_class'] for x in v}; return 'KNOWN_PV' if 'KNOWN_PV' in o else 'STRONG_PV_CANDIDATE' if 'STRONG_PV_CANDIDATE' in o else 'REPEAT_ONLY'
rows = []; elig = {}
for l, v in sorted(loci.items()):
    ce = any(x['caller_eligible'].startswith('YES') for x in v); p = plat.get(l)
    if not ce: st, why, uses = 'NOT_APPLICABLE_UNCALLABLE_OR_INELIGIBLE', 'not caller eligible (no uniquely anchorable member)', 'NONE'
    elif p and p['callability_class'] == 'PLATFORM_CONFOUNDED_CALLABILITY': st, why, uses = 'INELIGIBLE_PLATFORM_CONFOUNDED_CALLABILITY', 'callable-fraction difference of at least 0.5 between PRJNA360417 and PRJNA1103397 (PA02)', 'TECHNICAL_VALIDATION;WITHIN_COHORT_DESCRIPTION;PRESPECIFIED_SENSITIVITY_ANALYSIS'
    elif p and p['callability_class'] == 'UNCALLABLE_BOTH': st, why, uses = 'NOT_APPLICABLE_UNCALLABLE_OR_INELIGIBLE', 'callable fraction below 0.2 in both pilot studies', 'TECHNICAL_VALIDATION'
    elif p and p['callability_class'] == 'NA': st, why, uses = 'NOT_APPLICABLE_UNCALLABLE_OR_INELIGIBLE', 'callability not assessable in both pilot studies', 'TECHNICAL_VALIDATION'
    else: st, why, uses = 'PENDING_G7', 'not platform confounded; eligibility for disease analysis is decided in G7 (not automatic)', 'TECHNICAL_VALIDATION;WITHIN_COHORT_DESCRIPTION;PRESPECIFIED_SENSITIVITY_ANALYSIS;G7_ELIGIBILITY_REVIEW'
    elig[l] = st
    rows.append(dict(locus_id=l, evidence_class=cls(v), caller_eligible='YES' if ce else 'NO', callability_class=p['callability_class'] if p else 'NA', primary_association_eligibility=st, reason=why, allowed_uses=uses,
                     restoration_rule='documented disease-blind technical rule supported by G6/G7 evidence only (PP03); never via disease association'))
with open('metadata/g5/g5_locus_eligibility.tsv', 'w', newline='') as f:
    w = csv.DictWriter(f, list(rows[0]), delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(rows)
cols = [c for c in cat[0] if c != 'primary_association_eligibility'] + ['primary_association_eligibility']
for r in cat: r['primary_association_eligibility'] = elig[r['locus_id']]
with open('metadata/g5/g5_repeat_catalogue.tsv', 'w', newline='') as f:
    w = csv.DictWriter(f, cols, delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(cat)
print(collections.Counter(elig.values()))
