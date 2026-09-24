#!/usr/bin/env python3
"""G5 synthetic benchmark on the real catalogue (26695 + J99 members simulated on their own genomes; extraction and calling
use the FULL catalogue index so cross-member ambiguity is exercised). Stratified deterministic subset (no disease data).
usage: g5_benchmark_real.py RULES OUTPREFIX [NOISE_TSV]"""
import csv, gzip, os, sys, collections
sys.path.insert(0, 'scripts/g5')
from g5caller.common import Params, load_catalogue, load_params
from g5caller.benchmark import run_benchmark
from g5caller.noise import NoiseModel
rules, outdir = sys.argv[1], sys.argv[2]
mode = sys.argv[3] if len(sys.argv) > 3 else 'PLACEHOLDER'  # PLACEHOLDER | EMPIRICAL (stutter from PRJNA622860 invariant loci) | STRESS (2x empirical)
params = load_params(rules)
panel = 'data/reference/g5_panel/panel_all.fna.gz'  # all 8 panel genomes (built by g5_merge step)
cat = 'metadata/g5/g5_repeat_catalogue.tsv'
R = list(csv.DictReader(open(cat), delimiter='\t'))
def stratum(r):
    ml = int(r['motif_length']); base = r['motif_canonical']
    mc = 'MONO_AT' if base == 'A' else 'MONO_GC' if ml == 1 else f'MOTIF{ml}'
    return (mc, 'CDS' if r['context'] == 'CDS' else 'NONCDS')
strata = collections.defaultdict(list)
ORDER = {'26695': 0, 'J99': 1}
for r in sorted(R, key=lambda r: (ORDER.get(r['reference_short'], 2), r['member_id'])):
    # 26695 and J99 members first; other panel genomes only fill strata that would otherwise have fewer than PER members
    if r['caller_eligible'].startswith('YES'): strata[stratum(r)].append(r)
PER = 10
ids = set(); rows = []; used_loci = set()
for k in sorted(strata):
    v = strata[k]; main = [r for r in v if r['reference_short'] in ('26695', 'J99')]; rest = [r for r in v if r['reference_short'] not in ('26695', 'J99')]
    step = max(1, len(main) // PER); pick = main[::step][:PER]
    if len(pick) < PER: pick += rest[::max(1, len(rest) // (PER - len(pick)))][:PER - len(pick)]
    for r in pick:
        if r['locus_id'] not in used_loci: ids.add(r['member_id']); used_loci.add(r['locus_id'])  # one member per locus: orthologs of one locus cannot carry different alleles in one sample
for r in R:  # always include KNOWN_PV eligible members of the two curated references
    if r['evidence_class'] == 'KNOWN_PV' and r['reference_short'] in ('26695', 'J99') and r['caller_eligible'].startswith('YES') and r['locus_id'] not in used_loci: ids.add(r['member_id']); used_loci.add(r['locus_id'])
print('simulated members', len(ids), {k: len(v) for k, v in sorted(strata.items())}, flush=True)
from g5caller.simulate import StutterSpec
spec = None
if mode in ('EMPIRICAL', 'STRESS'):
    est = {r['param_id']: float(r['value']) for r in csv.DictReader(open('results/g5/work/noise_622860_estimated.tsv'), delimiter='\t')}
    f = 2.0 if mode == 'STRESS' else 1.0
    spec = StutterSpec({k: f * est['stutter_a_' + k] for k in ('MONO', 'DI', 'POLY')}, {k: f * est['stutter_b_' + k] for k in ('MONO', 'DI', 'POLY')}, est['stutter_up_ratio'], est['stutter_k2_ratio'])
noise = None
per = []
byref = collections.defaultdict(set)
for r in R:
    if r['member_id'] in ids: byref[r['reference_short']].add(r['member_id'])
for ref in sorted(byref):  # one genome per synthetic sample (a real sample carries one genome); extraction/calling use the FULL catalogue
    p_, s_ = run_benchmark(cat, panel, os.path.join(outdir, ref), layouts=('SE150', 'PE250'), depths=(10, 30, 100, 300), mixes=('100:0', '90:10', '70:30', '50:50'), stutter=True, params=params, noise=noise, only_ids=byref[ref], stutter_spec=spec)
    per += p_
    print('genome', ref, len(byref[ref]), 'members', flush=True)
info = {r['member_id']: r for r in R}
for r in per:
    i = info[r['member_id']]; r['stratum'] = '%s|%s' % stratum(i); r['evidence_class'] = i['evidence_class']; r['reference_short'] = i['reference_short']
    r['tract_bp'] = i['tract_length_bp']
with open(os.path.join(outdir, 'g5_synthetic_benchmark_per_locus.tsv'), 'w', newline='') as f:
    w = csv.DictWriter(f, list(per[0]), delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(per)
print('done', len(per))
