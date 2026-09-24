#!/usr/bin/env python3
"""G7 section 5: descriptive confounding audit of disease stage versus
study, country, platform, biopsy site, patient multiplicity and (where
sequence data allow) population structure. Purely descriptive counts; does
not test any phase-variable locus against disease and does not declare a
disease effect (CLAUDE.md rule 10, G7 issue section 5).

Source: metadata/g3/g3_candidate_biosamples.tsv (the full G3 harmonized
candidate biosample table, 1,884 rows, 793 with a Correa-stage disease
label across 44 BioProjects). This table already carries verification_status
per BioSample from G1/G2/G3; that status is preserved here and NOT upgraded.
Only Correa-labelled rows (NAG/AG/IM/DYS/GC) are tabulated; UNKNOWN/OTHER
rows are excluded from the contingency tables but counted in the summary.

Run: python3 scripts/g7/g7_confounding_audit.py
"""
import csv
import os
from collections import Counter, defaultdict

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SRC = os.path.join(REPO, "metadata", "g3", "g3_candidate_biosamples.tsv")
STRUCT_COORDS = os.path.join(REPO, "results", "g7", "g7_structure_coordinates.tsv")
OUT_DIR = os.path.join(REPO, "results", "g7")

STAGES = ["NAG", "AG", "IM", "DYS", "GC"]


def load():
    with open(SRC, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_contingency(path, rows, row_key, col_key=None):
    """row_key groups rows (e.g. bioproject); columns are the 5 Correa
    stages plus totals and verification status breakdown."""
    groups = defaultdict(lambda: Counter())
    verif = defaultdict(lambda: Counter())
    country = defaultdict(lambda: set())
    for r in rows:
        if r["disease_harmonized"] not in STAGES:
            continue
        k = r[row_key]
        groups[k][r["disease_harmonized"]] += 1
        verif[k][r["verification_status"]] += 1
        country[k].add(r["country_harmonized"])
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow([row_key] + STAGES + ["total_correa", "n_verification_statuses", "verification_status_counts", "n_countries", "countries"])
        for k in sorted(groups.keys()):
            c = groups[k]
            total = sum(c.values())
            vstr = ";".join(f"{s}:{n}" for s, n in sorted(verif[k].items()))
            w.writerow([k] + [c.get(s, 0) for s in STAGES] + [total, len(verif[k]), vstr, len(country[k]), ";".join(sorted(country[k]))])


def classify_contrast(counts_by_group, group_key_name):
    """Very simple, transparent identifiability rule, applied consistently:
    a stage-vs-group contrast is NOT_IDENTIFIABLE if fewer than 2 groups
    have >=1 sample in at least 2 different Correa stages with
    verification_status not REPOSITORY_ONLY/UNRESOLVED (i.e. at least
    SOURCE_VERIFIED_PARTIAL or better); EXPLORATORY_ONLY if such overlap
    exists but only in REPOSITORY_ONLY / low-confidence rows; otherwise the
    contrast is at least described as within-study/cross-study candidate
    (final defensibility judged qualitatively in the G7 report, not
    auto-declared DEFENSIBLE here)."""
    return None  # descriptive helper retained for documentation; classification done in report


def main():
    rows = load()
    correa = [r for r in rows if r["disease_harmonized"] in STAGES]

    write_contingency(os.path.join(OUT_DIR, "g7_stage_by_study.tsv"), rows, "bioproject")
    write_contingency(os.path.join(OUT_DIR, "g7_stage_by_country.tsv"), rows, "country_harmonized")

    # platform: platform_classes can be multi-valued (e.g. "ILLUMINA_SHORT")
    write_contingency(os.path.join(OUT_DIR, "g7_stage_by_platform.tsv"), rows, "platform_classes")

    # structure: only the 32-sample population-structure dataset has coordinates;
    # join by biosample id where possible (PRJNA521871 only)
    struct = {}
    if os.path.exists(STRUCT_COORDS):
        with open(STRUCT_COORDS, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                struct[r["sample_id"]] = r
    with open(os.path.join(OUT_DIR, "g7_stage_by_structure.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["biosample", "bioproject", "disease_harmonized", "PCo1", "PCo2", "PCo3", "structure_role"])
        n_joined = 0
        for r in correa:
            s = struct.get(r["biosample"])
            if s:
                n_joined += 1
                w.writerow([r["biosample"], r["bioproject"], r["disease_harmonized"], s["PCo1"], s["PCo2"], s["PCo3"], s["role"]])
        print("stage_by_structure: disease-labelled BioSamples with computed structure coordinates:", n_joined, "of", len(correa))

    # overall confounding matrix: one row per bioproject with disease x country x platform x patient-mapping summary
    by_study = defaultdict(list)
    for r in correa:
        by_study[r["bioproject"]].append(r)
    with open(os.path.join(OUT_DIR, "g7_confounding_matrix.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["bioproject", "n_correa_biosamples", "stage_counts", "countries", "platform_classes",
                     "verification_status_counts", "patient_mapping_status_counts", "n_distinct_patients_confirmed",
                     "public_reads_status_counts", "biopsy_site_counts", "has_structure_coordinates"])
        for study, srows in sorted(by_study.items()):
            stage_c = Counter(r["disease_harmonized"] for r in srows)
            countries = sorted(set(r["country_harmonized"] for r in srows))
            platforms = sorted(set(r["platform_classes"] for r in srows))
            verif_c = Counter(r["verification_status"] for r in srows)
            pm_c = Counter(r["patient_mapping_status"] for r in srows)
            confirmed_patients = len(set(r["patient_uid"] for r in srows if r["patient_mapping_status"] == "CONFIRMED" and r["patient_uid"] not in ("NA", "")))
            reads_c = Counter(r["public_reads_status"] for r in srows)
            site_c = Counter(r["biopsy_site_harmonized"] for r in srows)
            has_struct = "YES" if any(r["biosample"] in struct for r in srows) else "NO"
            w.writerow([
                study, len(srows),
                ";".join(f"{s}:{n}" for s, n in sorted(stage_c.items())),
                ";".join(countries),
                ";".join(platforms),
                ";".join(f"{s}:{n}" for s, n in sorted(verif_c.items())),
                ";".join(f"{s}:{n}" for s, n in sorted(pm_c.items())),
                confirmed_patients,
                ";".join(f"{s}:{n}" for s, n in sorted(reads_c.items())),
                ";".join(f"{s}:{n}" for s, n in sorted(site_c.items())),
                has_struct,
            ])

    print("total candidate biosamples:", len(rows))
    print("Correa-labelled candidate biosamples:", len(correa))
    print("distinct studies with Correa label:", len(by_study))


if __name__ == "__main__":
    main()
