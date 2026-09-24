#!/usr/bin/env python3
"""G7 section 1: freeze the primary technical locus universe before any
phase-state vs disease summary.

Source of truth: metadata/g6/g6_primary_eligibility_by_locus.tsv, produced
and committed at G6 final closure under the frozen G6 policies
(metadata/g6/g6_primary_eligibility_policy.tsv, g5_platform_policy.tsv,
g6_boundary_audit_rules.tsv). This script does not recompute or alter any
G5/G6 threshold; it only relabels the already-frozen per-locus eligibility
into the three G7 classes named in the G7 issue (PRIMARY_TECHNICAL,
PROVISIONAL_TECHNICAL, PRIMARY_INELIGIBLE) and writes a frozen ruleset.

Run: python3 scripts/g7/g7_build_locus_universe.py
"""
import csv
import hashlib
import os

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SRC = os.path.join(REPO, "metadata", "g6", "g6_primary_eligibility_by_locus.tsv")
OUT_UNIVERSE = os.path.join(REPO, "metadata", "g7", "g7_locus_analysis_universe.tsv")
OUT_RULES = os.path.join(REPO, "metadata", "g7", "g7_analysis_rules.tsv")
OUT_RULES_SHA = os.path.join(REPO, "metadata", "g7", "g7_analysis_rules.sha256")

RULES = [
    ("G7L01", "PRIMARY_TECHNICAL = loci with g6_primary_eligibility_by_locus.g7_technical_eligibility == ELIGIBLE. "
               "These are the 261 mono A/T loci that are G6 BOUNDARY_STABLE, not G5 PLATFORM_CONFOUNDED_CALLABILITY, "
               "caller-eligible/callable under frozen G5 rules, and motif class EMPIRICALLY_VALIDATED_EXACT."),
    ("G7L02", "PROVISIONAL_TECHNICAL = loci with g7_technical_eligibility == ELIGIBLE_FLAGGED_PROVISIONAL. "
               "The 39 remaining G6 technically eligible loci (18 trinucleotide, 16 mono G/C, 5 dinucleotide), "
               "motif class SUPPORTED_BUT_ALLELE_DIVERSITY_LIMITED. These proceed through G7 technical/confounding "
               "audit but stay separately flagged and cannot alone support a strong primary disease claim without "
               "later orthogonal validation if they become key candidates."),
    ("G7L03", "PRIMARY_INELIGIBLE = every other locus in the 1,405-locus G5/G6 catalogue: all non-BOUNDARY_STABLE "
               "loci (BOUNDARY_AMBIGUOUS, ORTHOLOGY_AMBIGUOUS, INSUFFICIENT_TO_RESOLVE), all G5 "
               "PLATFORM_CONFOUNDED_CALLABILITY loci, all 8 KNOWN_PV loci as a primary class, all MIXTURE states, "
               "and all loci whose motif-class interpretation is INSUFFICIENT_EMPIRICAL_VALIDATION or that are not "
               "caller-eligible/callable under the frozen G5 rules."),
    ("G7L04", "No locus is rescued into PRIMARY_TECHNICAL or PROVISIONAL_TECHNICAL based on disease data; this "
               "freeze happens before any phase-state vs disease summary is computed (CLAUDE.md rule 10)."),
    ("G7L05", "This freeze does not modify any G5 caller threshold or any G6 eligibility/boundary-audit rule "
               "(g5_caller_rules.tsv, g6_boundary_audit_rules.tsv, g6_primary_eligibility_policy.tsv unchanged)."),
    ("G7L06", "Any future catalogue rescue of PRIMARY_INELIGIBLE loci must be separately pre-specified, "
               "disease-blind, and secondary/sensitivity-only unless explicitly re-authorized before disease-"
               "association inspection (supervisor decision recorded at G7 authorization)."),
]


def classify(row):
    elig = row["g7_technical_eligibility"]
    if elig == "ELIGIBLE":
        return "PRIMARY_TECHNICAL"
    if elig == "ELIGIBLE_FLAGGED_PROVISIONAL":
        return "PROVISIONAL_TECHNICAL"
    return "PRIMARY_INELIGIBLE"


def main():
    with open(SRC, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    out_fields = [
        "locus_id", "motif_class", "motif_class_interpretation", "known_pv",
        "g6_audit_status", "g5_platform_confounded_or_uncallable",
        "g7_locus_class", "g7_source_eligibility_value", "reasons",
    ]
    universe_rows = []
    for r in rows:
        g7_class = classify(r)
        reasons = r["reasons"]
        platform_or_uncallable = "YES" if (
            "PLATFORM_CONFOUNDED_CALLABILITY" in reasons
            or "NOT_CALLABLE_OR_NOT_CALLER_ELIGIBLE_UNDER_G5" in reasons
        ) else "NO"
        universe_rows.append({
            "locus_id": r["locus_id"],
            "motif_class": r["motif_class"],
            "motif_class_interpretation": r["motif_class_interpretation"],
            "known_pv": r["known_pv"],
            "g6_audit_status": r["audit_status"],
            "g5_platform_confounded_or_uncallable": platform_or_uncallable,
            "g7_locus_class": g7_class,
            "g7_source_eligibility_value": r["g7_technical_eligibility"],
            "reasons": reasons if reasons else "NA",
        })

    os.makedirs(os.path.dirname(OUT_UNIVERSE), exist_ok=True)
    with open(OUT_UNIVERSE, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=out_fields, delimiter="\t")
        w.writeheader()
        for r in universe_rows:
            w.writerow(r)

    with open(OUT_RULES, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["rule_id", "rule"])
        for rid, text in RULES:
            w.writerow([rid, text])

    sha = hashlib.sha256(open(OUT_RULES, "rb").read()).hexdigest()
    with open(OUT_RULES_SHA, "w") as fh:
        fh.write(f"{sha}  g7_analysis_rules.tsv\n")

    from collections import Counter
    counts = Counter(r["g7_locus_class"] for r in universe_rows)
    motif_by_class = {}
    for cls in ("PRIMARY_TECHNICAL", "PROVISIONAL_TECHNICAL"):
        motif_by_class[cls] = Counter(r["motif_class"] for r in universe_rows if r["g7_locus_class"] == cls)

    print("total loci:", len(universe_rows))
    print("counts:", dict(counts))
    print("PRIMARY_TECHNICAL by motif:", dict(motif_by_class["PRIMARY_TECHNICAL"]))
    print("PROVISIONAL_TECHNICAL by motif:", dict(motif_by_class["PROVISIONAL_TECHNICAL"]))
    print("rules sha256:", sha)


if __name__ == "__main__":
    main()
