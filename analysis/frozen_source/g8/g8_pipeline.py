#!/usr/bin/env python3
"""G8 pilot disease association (within-study only). Pure-Python (no numpy/scipy/pandas installed
in this environment, consistent with G7/G7-follow-up). Implements metadata/g8/g8_analysis_rules.tsv
(frozen by scripts/g8/g8_build_rules.py before this script ever reads a disease label). All inputs
are already-committed, already disease-blind-generated caller outputs from prior gates:
- results/g5/g5_phase1_calls.tsv (PRJNA360417, PRJNA1103397; full depth, G5 caller)
- results/g7_followup/g7f_locus_level_calls.tsv (PRJNA678459; 500,000-fragment bounded G5 caller,
  fixed exploratory input per Issue #9 / issue #12 instruction; NOT re-run here)
No raw read is opened by this script. No network access. No threshold in g8_analysis_rules.tsv is
changed by this script.
"""
import csv
import hashlib
import itertools
import json
import os
import random
import statistics

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.chdir(REPO)

SEED = 88022026
N_BOOT = 2000

RULES_PATH = 'metadata/g8/g8_analysis_rules.tsv'
RULES_SHA_PATH = RULES_PATH + '.sha256'


def verify_rules_frozen():
    digest = hashlib.sha256(open(RULES_PATH, 'rb').read()).hexdigest()
    expected = open(RULES_SHA_PATH).read().split()[0]
    assert digest == expected, 'g8_analysis_rules.tsv has drifted from its frozen hash'
    return digest


def load_locus_universe():
    primary, provisional = set(), set()
    with open('metadata/g7/g7_locus_analysis_universe.tsv') as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            if row['g7_locus_class'] == 'PRIMARY_TECHNICAL':
                primary.add(row['locus_id'])
            elif row['g7_locus_class'] == 'PROVISIONAL_TECHNICAL':
                provisional.add(row['locus_id'])
    assert len(primary) == 261, f'PRIMARY_TECHNICAL locus count drifted: {len(primary)}'
    assert len(provisional) == 39, f'PROVISIONAL_TECHNICAL locus count drifted: {len(provisional)}'
    return primary, provisional


def load_pilot_manifest():
    rows = list(csv.DictReader(open('metadata/g3/g3_disease_pilot_manifest.tsv'), delimiter='\t'))
    by_biosample = {r['biosample']: r for r in rows}
    return by_biosample


PRJNA678459_LABELS = {
    'SAMN16801845': 'AG', 'SAMN16801850': 'GC', 'SAMN16801846': 'AG', 'SAMN16801843': 'GC',
    'SAMN16801852': 'AG', 'SAMN16801851': 'GC', 'SAMN16801849': 'GC', 'SAMN16801848': 'GC',
    'SAMN16801844': 'AG', 'SAMN16801847': 'AG',
}  # metadata/harmonized/biosample_metadata.tsv, BioSample.host_disease, EXPLICIT, D02/D05 (verified this gate)


def load_g5_calls(target_biosamples, target_loci):
    calls = {}
    with open('results/g5/g5_phase1_calls.tsv') as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            bs, loc = row['biosample'], row['locus_id']
            if bs not in target_biosamples or loc not in target_loci:
                continue
            dom = row['dominant_bp']
            calls[(bs, loc)] = dict(
                callable=(row['callable'] == 'YES'),
                dominant_bp=(int(dom) if dom not in ('', 'NA', '.') else None),
            )
    return calls


def load_g7f_calls(target_biosamples, target_loci):
    calls = {}
    with open('results/g7_followup/g7f_locus_level_calls.tsv') as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            bs, loc = row['biosample'], row['locus_id']
            if bs not in target_biosamples or loc not in target_loci:
                continue
            dom = row['dominant_bp']
            calls[(bs, loc)] = dict(
                callable=(row['callable'] == 'YES'),
                dominant_bp=(int(dom) if dom not in ('', 'NA', '.') else None),
            )
    return calls


def screen_loci(calls, samples, loci):
    """R09: disease-blind screen. Returns dict locus_id -> dict(pass, callability, n_distinct, minor_count, reason)."""
    out = {}
    n = len(samples)
    for loc in loci:
        vals = []
        n_callable = 0
        for s in samples:
            c = calls.get((s, loc))
            if c and c['callable'] and c['dominant_bp'] is not None:
                n_callable += 1
                vals.append(c['dominant_bp'])
        callability = n_callable / n
        counts = {}
        for v in vals:
            counts[v] = counts.get(v, 0) + 1
        n_distinct = len(counts)
        minor_count = sorted(counts.values(), reverse=True)[1] if n_distinct >= 2 else 0
        ok = callability >= 0.80 and n_distinct >= 2 and minor_count >= 2
        reasons = []
        if callability < 0.80:
            reasons.append(f'CALLABILITY_{callability:.3f}_LT_0.80')
        if n_distinct < 2:
            reasons.append('MONOMORPHIC_OR_UNCALLABLE')
        if n_distinct >= 2 and minor_count < 2:
            reasons.append(f'MINOR_ALLELE_COUNT_{minor_count}_LT_2')
        out[loc] = dict(pass_screen=ok, n_callable=n_callable, callability=round(callability, 4),
                         n_distinct=n_distinct, minor_count=minor_count,
                         reason=(';'.join(reasons) if reasons else 'PASS'))
    return out


MIN_COVERAGE_PER_PAIR = 10  # R11 pre-specified floor; frozen before any distance value was inspected


def pairwise_distance_matrix(calls, samples, loci):
    """R11: fraction of jointly callable loci with different dominant_bp.
    Returns (D_all, coverage) where D_all[pair] is always defined whenever joint>=1 (diagnostic use);
    callers must separately check coverage[pair] >= MIN_COVERAGE_PER_PAIR before treating a distance
    as meeting the pre-specified primary-analysis floor."""
    D = {}
    coverage = {}
    for a, b in itertools.combinations(samples, 2):
        joint = 0
        diff = 0
        for loc in loci:
            ca = calls.get((a, loc))
            cb = calls.get((b, loc))
            if ca and cb and ca['callable'] and cb['callable'] and ca['dominant_bp'] is not None and cb['dominant_bp'] is not None:
                joint += 1
                if ca['dominant_bp'] != cb['dominant_bp']:
                    diff += 1
        key = frozenset((a, b))
        coverage[key] = joint
        D[key] = (diff / joint) if joint >= 1 else None
    return D, coverage


def permanova_stat(D, ids, groups):
    n = len(ids)
    total_sq = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            d = D[frozenset((ids[i], ids[j]))]
            total_sq += d * d
    SST = total_sq / n
    SSW = 0.0
    for members in groups.values():
        ng = len(members)
        s = 0.0
        for i in range(ng):
            for j in range(i + 1, ng):
                d = D[frozenset((members[i], members[j]))]
                s += d * d
        SSW += s / ng
    k = len(groups)
    SSA = SST - SSW
    F = (SSA / (k - 1)) / (SSW / (n - k)) if (SSW > 0 and n > k) else float('nan')
    R2 = SSA / SST if SST > 0 else float('nan')
    return F, R2


def exhaustive_permanova_test(D, group_a_ids, group_b_ids):
    ids = group_a_ids + group_b_ids
    n1 = len(group_a_ids)
    n = len(ids)
    obs_groups = {'A': list(group_a_ids), 'B': list(group_b_ids)}
    F_obs, R2_obs = permanova_stat(D, ids, obs_groups)
    total = 0
    ge = 0
    for combo in itertools.combinations(ids, n1):
        setA = list(combo)
        setB = [x for x in ids if x not in set(combo)]
        F_p, _ = permanova_stat(D, ids, {'A': setA, 'B': setB})
        total += 1
        if F_p >= F_obs - 1e-12:
            ge += 1
    p_exact = ge / total
    return dict(pseudo_F=F_obs, R2=R2_obs, exact_p=p_exact, n_permutations=total)


def bootstrap_r2_ci(D, group_a_ids, group_b_ids, seed=SEED, n_boot=N_BOOT):
    rng = random.Random(seed)
    r2s = []
    for _ in range(n_boot):
        drawA = [rng.choice(group_a_ids) for _ in range(len(group_a_ids))]
        drawB = [rng.choice(group_b_ids) for _ in range(len(group_b_ids))]
        # build pseudo-sample positional labels to allow repeats; distance between two positions
        # drawing the same original id is 0 (identical underlying call vector); otherwise the
        # original pairwise distance is reused (case resampling bootstrap on the distance matrix).
        posA = [f'A{i}' for i in range(len(drawA))]
        posB = [f'B{i}' for i in range(len(drawB))]
        orig = dict(zip(posA, drawA))
        orig.update(dict(zip(posB, drawB)))
        allpos = posA + posB

        class LazyD(dict):
            def __getitem__(self, key):
                p, q = tuple(key)
                op, oq = orig[p], orig[q]
                if op == oq:
                    return 0.0
                return D[frozenset((op, oq))]
        r2 = permanova_stat(LazyD(), allpos, {'A': posA, 'B': posB})[1]
        if r2 == r2:  # skip nan
            r2s.append(r2)
    r2s.sort()
    lo = r2s[int(0.025 * len(r2s))]
    hi = r2s[int(0.975 * len(r2s)) - 1]
    return round(lo, 4), round(hi, 4), len(r2s)


def cliffs_delta(a_vals, b_vals):
    gt = lt = 0
    for x in a_vals:
        for y in b_vals:
            if y > x:
                gt += 1
            elif y < x:
                lt += 1
    n = len(a_vals) * len(b_vals)
    return (gt - lt) / n if n else float('nan')


def exhaustive_mean_diff_test(a_vals, b_vals):
    ids = list(range(len(a_vals) + len(b_vals)))
    vals = a_vals + b_vals
    n1 = len(a_vals)
    obs = statistics.mean(b_vals) - statistics.mean(a_vals)
    total = 0
    ge = 0
    for combo in itertools.combinations(ids, n1):
        setA = [vals[i] for i in combo]
        setB = [vals[i] for i in ids if i not in set(combo)]
        stat = abs(statistics.mean(setB) - statistics.mean(setA))
        total += 1
        if stat >= abs(obs) - 1e-12:
            ge += 1
    return obs, ge / total, total


def bootstrap_meandiff_ci(a_vals, b_vals, seed=SEED, n_boot=N_BOOT):
    rng = random.Random(seed)
    diffs = []
    for _ in range(n_boot):
        ra = [rng.choice(a_vals) for _ in a_vals]
        rb = [rng.choice(b_vals) for _ in b_vals]
        diffs.append(statistics.mean(rb) - statistics.mean(ra))
    diffs.sort()
    lo = diffs[int(0.025 * len(diffs))]
    hi = diffs[int(0.975 * len(diffs)) - 1]
    return round(lo, 4), round(hi, 4)


def bh_qvalues(pvals):
    m = len(pvals)
    order = sorted(range(m), key=lambda i: pvals[i])
    q = [0.0] * m
    prev = 1.0
    for rank, idx in enumerate(reversed(order), start=1):
        i = m - rank + 1
        p = pvals[idx]
        val = min(prev, p * m / i)
        prev = val
        q[idx] = val
    return q


def locus_level_table(calls, group_a, group_b, loci, screen):
    rows = []
    pvals = []
    for loc in loci:
        sc = screen[loc]
        if not sc['pass_screen']:
            continue
        a_vals = [calls[(s, loc)]['dominant_bp'] for s in group_a
                  if calls.get((s, loc)) and calls[(s, loc)]['callable'] and calls[(s, loc)]['dominant_bp'] is not None]
        b_vals = [calls[(s, loc)]['dominant_bp'] for s in group_b
                  if calls.get((s, loc)) and calls[(s, loc)]['callable'] and calls[(s, loc)]['dominant_bp'] is not None]
        if len(a_vals) < 2 or len(b_vals) < 2:
            continue
        obs_diff, p_exact, n_perm = exhaustive_mean_diff_test(a_vals, b_vals)
        delta = cliffs_delta(a_vals, b_vals)
        lo, hi = bootstrap_meandiff_ci(a_vals, b_vals)
        rows.append(dict(
            locus_id=loc, n_a=len(a_vals), n_b=len(b_vals),
            median_a=statistics.median(a_vals), median_b=statistics.median(b_vals),
            mean_diff_b_minus_a=round(obs_diff, 4), cliffs_delta=round(delta, 4),
            boot_ci_lo=lo, boot_ci_hi=hi, exact_p=round(p_exact, 6), n_permutations=n_perm,
        ))
        pvals.append(p_exact)
    if rows:
        qs = bh_qvalues(pvals)
        for r, q in zip(rows, qs):
            r['bh_q'] = round(q, 6)
    return rows


def cohort_level_endpoint(D_all, coverage, group_a, group_b, n_loci_available):
    """Computes the R11/R12 cohort-level endpoint two ways:
    - STRICT: only valid if every sample pair meets the pre-specified MIN_COVERAGE_PER_PAIR floor
      (R11); otherwise NOT_COMPUTABLE, disclosed rather than silently downgraded.
    - DIAGNOSTIC: uses whatever coverage exists for every pair (no floor), explicitly labeled
      low-coverage/diagnostic-only when n_loci_available < MIN_COVERAGE_PER_PAIR; never used alone to
      claim a positive or negative result."""
    ids = group_a + group_b
    all_pairs = [frozenset(p) for p in itertools.combinations(ids, 2)]
    min_cov = min(coverage[p] for p in all_pairs)
    strict_ok = min_cov >= MIN_COVERAGE_PER_PAIR
    result = dict(n_loci_available=n_loci_available, min_pairwise_coverage=min_cov,
                  max_pairwise_coverage=max(coverage[p] for p in all_pairs),
                  strict_computable=strict_ok)
    diag_D = {p: D_all[p] for p in all_pairs if D_all[p] is not None}
    if len(diag_D) == len(all_pairs):
        perm = exhaustive_permanova_test(diag_D, group_a, group_b)
        boot = bootstrap_r2_ci(diag_D, group_a, group_b)
        result.update(diagnostic_pseudo_F=round(perm['pseudo_F'], 4), diagnostic_R2=round(perm['R2'], 4),
                      diagnostic_exact_p=round(perm['exact_p'], 6), diagnostic_n_permutations=perm['n_permutations'],
                      diagnostic_boot_r2_ci_lo=boot[0], diagnostic_boot_r2_ci_hi=boot[1],
                      diagnostic_label=('BELOW_PRESPECIFIED_COVERAGE_FLOOR_DIAGNOSTIC_ONLY' if not strict_ok else 'MEETS_COVERAGE_FLOOR'))
        if strict_ok:
            result.update(strict_pseudo_F=result['diagnostic_pseudo_F'], strict_R2=result['diagnostic_R2'],
                          strict_exact_p=result['diagnostic_exact_p'], strict_n_permutations=result['diagnostic_n_permutations'],
                          strict_boot_r2_ci_lo=result['diagnostic_boot_r2_ci_lo'], strict_boot_r2_ci_hi=result['diagnostic_boot_r2_ci_hi'])
        else:
            result.update(strict_pseudo_F='NOT_COMPUTABLE', strict_R2='NOT_COMPUTABLE', strict_exact_p='NOT_COMPUTABLE',
                          strict_n_permutations='NOT_COMPUTABLE', strict_boot_r2_ci_lo='NOT_COMPUTABLE', strict_boot_r2_ci_hi='NOT_COMPUTABLE')
    else:
        result.update(strict_pseudo_F='NOT_COMPUTABLE_MISSING_PAIR', diagnostic_label='NOT_COMPUTABLE_MISSING_PAIR')
    return result


def write_tsv(path, rows, fields):
    with open(path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fields, delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def main():
    rules_sha = verify_rules_frozen()
    primary, provisional = load_locus_universe()
    manifest = load_pilot_manifest()

    p360 = [b for b, r in manifest.items() if r['cohort_study'] == 'PRJNA360417']
    p1103 = [b for b, r in manifest.items() if r['cohort_study'] == 'PRJNA1103397']
    p678459 = list(PRJNA678459_LABELS.keys())

    labels_360 = {b: manifest[b]['disease_stage_original_g2_label'] for b in p360}
    labels_1103 = {b: manifest[b]['disease_stage_original_g2_label'] for b in p1103}

    calls_360 = load_g5_calls(set(p360), primary | provisional)
    calls_1103 = load_g5_calls(set(p1103), primary | provisional)
    calls_678459 = load_g7f_calls(set(p678459), primary | provisional)

    summary = dict(rules_sha256=rules_sha, seed=SEED, n_bootstrap=N_BOOT)

    # ---------- PRJNA360417 NAG vs IM ----------
    nag = sorted([b for b, d in labels_360.items() if d == 'NAG'])
    im = sorted([b for b, d in labels_360.items() if d == 'IM'])
    ag360 = sorted([b for b, d in labels_360.items() if d == 'AG'])
    dys360 = sorted([b for b, d in labels_360.items() if d == 'DYS'])
    assert len(nag) == 5 and len(im) == 6 and len(ag360) == 1 and len(dys360) == 1

    contrast_samples = nag + im
    screen_360_primary = screen_loci(calls_360, contrast_samples, primary)
    screen_360_provisional = screen_loci(calls_360, contrast_samples, provisional)
    n_pass_primary_360 = sum(1 for v in screen_360_primary.values() if v['pass_screen'])
    n_pass_provisional_360 = sum(1 for v in screen_360_provisional.values() if v['pass_screen'])

    passed_primary_360 = [l for l, v in screen_360_primary.items() if v['pass_screen']]
    D_360, cov_360 = pairwise_distance_matrix(calls_360, contrast_samples, passed_primary_360)
    cohort_360 = cohort_level_endpoint(D_360, cov_360, nag, im, len(passed_primary_360))

    locus_360_primary = locus_level_table(calls_360, nag, im, passed_primary_360, screen_360_primary)
    passed_prov_360 = [l for l, v in screen_360_provisional.items() if v['pass_screen']]
    locus_360_provisional = locus_level_table(calls_360, nag, im, passed_prov_360, screen_360_provisional)

    # AG/DYS descriptive only (n=1 each): report their dominant_bp at the passed PRIMARY loci vs NAG/IM medians
    descriptive_360 = []
    for loc in passed_primary_360:
        row = dict(locus_id=loc)
        for grp_name, grp in [('NAG', nag), ('IM', im), ('AG', ag360), ('DYS', dys360)]:
            vals = [calls_360[(s, loc)]['dominant_bp'] for s in grp
                    if calls_360.get((s, loc)) and calls_360[(s, loc)]['callable'] and calls_360[(s, loc)]['dominant_bp'] is not None]
            row[f'{grp_name}_values'] = ';'.join(str(v) for v in vals) if vals else 'NA'
        descriptive_360.append(row)

    # ---------- PRJNA678459 AG vs GC ----------
    ag = sorted([b for b, d in PRJNA678459_LABELS.items() if d == 'AG'])
    gc = sorted([b for b, d in PRJNA678459_LABELS.items() if d == 'GC'])
    assert len(ag) == 5 and len(gc) == 5

    screen_678_primary = screen_loci(calls_678459, ag + gc, primary)
    screen_678_provisional = screen_loci(calls_678459, ag + gc, provisional)
    n_pass_primary_678 = sum(1 for v in screen_678_primary.values() if v['pass_screen'])
    n_pass_provisional_678 = sum(1 for v in screen_678_provisional.values() if v['pass_screen'])

    passed_primary_678 = [l for l, v in screen_678_primary.items() if v['pass_screen']]
    D_678, cov_678 = pairwise_distance_matrix(calls_678459, ag + gc, passed_primary_678)
    cohort_678 = cohort_level_endpoint(D_678, cov_678, ag, gc, len(passed_primary_678))

    locus_678_primary = locus_level_table(calls_678459, ag, gc, passed_primary_678, screen_678_primary)
    passed_prov_678 = [l for l, v in screen_678_provisional.items() if v['pass_screen']]
    locus_678_provisional = locus_level_table(calls_678459, ag, gc, passed_prov_678, screen_678_provisional)

    # ---------- PRJNA1103397 GC-only descriptive ----------
    gc1103 = sorted(p1103)
    assert len(gc1103) == 6 and all(labels_1103[b] == 'GC' for b in gc1103)
    screen_1103_primary = screen_loci(calls_1103, gc1103, primary)
    n_pass_primary_1103 = sum(1 for v in screen_1103_primary.values() if v['pass_screen'])
    passed_primary_1103 = [l for l, v in screen_1103_primary.items() if v['pass_screen']]
    descriptive_1103 = []
    for loc in passed_primary_1103:
        vals = [calls_1103[(s, loc)]['dominant_bp'] for s in gc1103
                if calls_1103.get((s, loc)) and calls_1103[(s, loc)]['callable'] and calls_1103[(s, loc)]['dominant_bp'] is not None]
        if len(vals) < 2:
            continue
        descriptive_1103.append(dict(locus_id=loc, n=len(vals), median=statistics.median(vals),
                                      distinct_values=len(set(vals)), values=';'.join(str(v) for v in vals)))

    # ---------- write outputs ----------
    def screen_rows(screen_dict, cohort, contrast, locus_class):
        return [dict(cohort=cohort, contrast=contrast, locus_class=locus_class, locus_id=l, **v) for l, v in screen_dict.items()]

    all_screen_rows = (
        screen_rows(screen_360_primary, 'PRJNA360417', 'NAG_vs_IM', 'PRIMARY_TECHNICAL') +
        screen_rows(screen_360_provisional, 'PRJNA360417', 'NAG_vs_IM', 'PROVISIONAL_TECHNICAL') +
        screen_rows(screen_678_primary, 'PRJNA678459', 'AG_vs_GC', 'PRIMARY_TECHNICAL') +
        screen_rows(screen_678_provisional, 'PRJNA678459', 'AG_vs_GC', 'PROVISIONAL_TECHNICAL') +
        screen_rows(screen_1103_primary, 'PRJNA1103397', 'GC_only_descriptive', 'PRIMARY_TECHNICAL')
    )
    write_tsv('results/g8/g8_locus_screen.tsv', all_screen_rows,
              ['cohort', 'contrast', 'locus_class', 'locus_id', 'pass_screen', 'n_callable', 'callability', 'n_distinct', 'minor_count', 'reason'])

    for rows, name in [(locus_360_primary, 'g8_locus_level_prjna360417_nag_vs_im_primary.tsv'),
                       (locus_360_provisional, 'g8_locus_level_prjna360417_nag_vs_im_provisional_sensitivity.tsv'),
                       (locus_678_primary, 'g8_locus_level_prjna678459_ag_vs_gc_primary.tsv'),
                       (locus_678_provisional, 'g8_locus_level_prjna678459_ag_vs_gc_provisional_sensitivity.tsv')]:
        if rows:
            write_tsv(f'results/g8/{name}', rows,
                      ['locus_id', 'n_a', 'n_b', 'median_a', 'median_b', 'mean_diff_b_minus_a',
                       'cliffs_delta', 'boot_ci_lo', 'boot_ci_hi', 'exact_p', 'n_permutations', 'bh_q'])
        else:
            open(f'results/g8/{name}', 'w').write('NO_LOCI_PASSED_SCREEN\n')

    write_tsv('results/g8/g8_descriptive_prjna360417_ag_dys_context.tsv', descriptive_360,
              ['locus_id', 'NAG_values', 'IM_values', 'AG_values', 'DYS_values'])
    write_tsv('results/g8/g8_descriptive_prjna1103397_gc_only.tsv', descriptive_1103,
              ['locus_id', 'n', 'median', 'distinct_values', 'values'])

    cohort_summary = [
        dict(cohort='PRJNA360417', contrast='NAG_vs_IM(5v6)', n_primary_passed_screen=n_pass_primary_360,
             n_provisional_passed_screen=n_pass_provisional_360, **cohort_360),
        dict(cohort='PRJNA678459', contrast='AG_vs_GC(5v5)_EXPLORATORY', n_primary_passed_screen=n_pass_primary_678,
             n_provisional_passed_screen=n_pass_provisional_678, **cohort_678),
        dict(cohort='PRJNA1103397', contrast='GC_only_descriptive_NO_TEST', n_primary_passed_screen=n_pass_primary_1103,
             n_provisional_passed_screen='NA', n_loci_available='NA', min_pairwise_coverage='NA',
             max_pairwise_coverage='NA', strict_computable='NA', strict_pseudo_F='NA', strict_R2='NA',
             strict_exact_p='NA', strict_n_permutations='NA', strict_boot_r2_ci_lo='NA', strict_boot_r2_ci_hi='NA',
             diagnostic_pseudo_F='NA', diagnostic_R2='NA', diagnostic_exact_p='NA', diagnostic_n_permutations='NA',
             diagnostic_boot_r2_ci_lo='NA', diagnostic_boot_r2_ci_hi='NA', diagnostic_label='NO_CONTRAST'),
    ]
    fields = ['cohort', 'contrast', 'n_primary_passed_screen', 'n_provisional_passed_screen', 'n_loci_available',
              'min_pairwise_coverage', 'max_pairwise_coverage', 'strict_computable',
              'strict_pseudo_F', 'strict_R2', 'strict_exact_p', 'strict_n_permutations',
              'strict_boot_r2_ci_lo', 'strict_boot_r2_ci_hi',
              'diagnostic_pseudo_F', 'diagnostic_R2', 'diagnostic_exact_p', 'diagnostic_n_permutations',
              'diagnostic_boot_r2_ci_lo', 'diagnostic_boot_r2_ci_hi', 'diagnostic_label']
    for row in cohort_summary:
        for f in fields:
            row.setdefault(f, 'NA')
    write_tsv('results/g8/g8_cohort_level_summary.tsv', cohort_summary, fields)

    summary.update(dict(
        prjna360417=dict(n_nag=5, n_im=6, n_ag_descriptive=1, n_dys_descriptive=1,
                          n_primary_screened_passed=n_pass_primary_360, n_provisional_screened_passed=n_pass_provisional_360,
                          n_locus_level_primary_tested=len(locus_360_primary), n_locus_level_provisional_tested=len(locus_360_provisional),
                          cohort_level=cohort_summary[0]),
        prjna678459=dict(n_ag=5, n_gc=5, exploratory_lower_provenance=True, caller_input='500000_fragments_per_sample_disease_blind_bounded_ISSUE9_output_NOT_RERUN',
                          n_primary_screened_passed=n_pass_primary_678, n_provisional_screened_passed=n_pass_provisional_678,
                          n_locus_level_primary_tested=len(locus_678_primary), n_locus_level_provisional_tested=len(locus_678_provisional),
                          cohort_level=cohort_summary[1]),
        prjna1103397=dict(n_gc=6, role='DESCRIPTIVE_CONTEXT_ONLY_NO_CONTRAST', n_primary_screened_passed=n_pass_primary_1103,
                           n_locus_descriptive_reported=len(descriptive_1103)),
    ))
    json.dump(summary, open('results/g8/g8_summary.json', 'w'), indent=2, default=str)
    print(json.dumps(summary, indent=2, default=str))


if __name__ == '__main__':
    main()
