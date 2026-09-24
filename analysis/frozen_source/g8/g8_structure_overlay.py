#!/usr/bin/env python3
"""G8 item 6: population-structure/confounding check, disease-label overlay only (structure itself
was computed disease-blind in prior gates; R19/R20). No new genome-wide distance is computed here.

PRJNA678459: extracts the 10x10 AG/GC submatrix from the already-frozen, disease-blind
results/g7_followup/g7f_structure_distance_matrix.tsv (read-based k-mer/MinHash sketches, bounded to
30,000 read pairs/sample, computed in the G7 follow-up gate before any G8 work existed) and overlays
the AG/GC labels descriptively: does the disease grouping look separated on this structure metric?
Also runs the same PERMANOVA machinery used for the repeat-allele-state distance endpoint on this structure distance,
purely as a descriptive comparison -- NOT a corrected/adjusted disease test and NOT a PV result.

PRJNA360417: NOT_PERFORMED. No existing disease-blind genome-wide/read-based structure result exists
for this study in the repository (checked: results/g7, results/g7_followup contain no PRJNA360417
entries), and the raw reads used for its G4/G5 Phase-1 processing are no longer present under
data/raw (checked this gate: only data/raw/g7_followup/PRJNA678459 exists). Recomputing from reads
would require RAW_READ_DOWNLOAD, which is not authorized in this gate (issue #12 lists only
GATE_TRANSITION G8). This gap is disclosed, not approximated or inferred.
"""
import csv
import itertools
import os
import sys

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.chdir(REPO)
sys.path.insert(0, 'scripts/g8')
from g8_pipeline import permanova_stat, exhaustive_permanova_test, PRJNA678459_LABELS  # noqa: E402


def main():
    rows = list(csv.reader(open('results/g7_followup/g7f_structure_distance_matrix.tsv'), delimiter='\t'))
    header = rows[0][1:]
    idx = {sid: i for i, sid in enumerate(header)}
    D = {}
    for r in rows[1:]:
        sid = r[0]
        for j, sid2 in enumerate(header):
            D[frozenset((sid, sid2))] = float(r[1 + j]) if sid != sid2 else 0.0

    ag = sorted([b for b, d in PRJNA678459_LABELS.items() if d == 'AG'])
    gc = sorted([b for b, d in PRJNA678459_LABELS.items() if d == 'GC'])
    ids = ag + gc

    sub = [['sample_id', 'disease_label'] + ids]
    for s in ids:
        lab = PRJNA678459_LABELS[s]
        sub.append([s, lab] + [f'{D[frozenset((s, t))]:.6f}' for t in ids])
    with open('results/g8/g8_structure_overlay_prjna678459_submatrix.tsv', 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n')
        w.writerows(sub)

    perm = exhaustive_permanova_test(D, ag, gc)

    nn_rows = list(csv.DictReader(open('results/g7_followup/g7f_structure_nearest_neighbors.tsv'), delimiter='\t'))
    nn_by_id = {r['biosample']: r['nearest_genome_1'] for r in nn_rows}
    all_same_nearest = len(set(nn_by_id.values())) == 1

    with open('results/g8/g8_structure_overlay_summary.tsv', 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n')
        w.writerow(['cohort', 'method', 'n_ag', 'n_gc', 'pseudo_F', 'R2', 'exact_p', 'n_permutations',
                    'all_10_samples_nearest_to_same_reference_genome', 'interpretation'])
        w.writerow(['PRJNA678459', 'SCREENING_STRUCTURE_ONLY_read_based_kmer_minhash_bounded_30000pairs_disease_blind_G7followup',
                    5, 5, round(perm['pseudo_F'], 4), round(perm['R2'], 4), round(perm['exact_p'], 6), perm['n_permutations'],
                    all_same_nearest,
                    'Descriptive only; no ancestry adjustment or definitive confounding correction claimed (R20). '
                    'AG/GC labels do not show obvious/perfect separation on this disease-blind screening-structure metric.'
                    if perm['exact_p'] > 0.05 else
                    'Descriptive only; AG/GC labels show a nominal association on this disease-blind screening-structure '
                    'metric (p<=0.05) -- this would be a confounding warning, not evidence of a PV-disease effect, and '
                    'is reported regardless of direction (R20, R22).'])
        w.writerow(['PRJNA360417', 'NOT_PERFORMED', 5, 6, 'NA', 'NA', 'NA', 'NA', 'NA',
                    'No existing disease-blind structure result for this study and its raw reads are no longer '
                    'local; RAW_READ_DOWNLOAD not authorized in this gate (R19). Disclosed gap, not approximated.'])

    print('PRJNA678459 structure-overlay PERMANOVA:', perm)
    print('all_same_nearest_reference_genome:', all_same_nearest)


if __name__ == '__main__':
    main()
