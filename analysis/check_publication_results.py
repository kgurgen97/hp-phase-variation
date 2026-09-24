#!/usr/bin/env python3
"""Fail if publication-facing tracked results drift from the frozen Project A record.

This script performs consistency checks only. It does not run the caller, download data,
change thresholds, or recompute biological associations.
"""
from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def read_tsv(path: str):
    with (ROOT / path).open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def require(condition: bool, message: str):
    if not condition:
        raise AssertionError(message)


def main():
    synthetic = json.loads((ROOT / "results/validation/synthetic_overall.json").read_text())
    s = synthetic["all"]
    require(s["callable"] == 1780, f"synthetic callable drift: {s['callable']}")
    require(abs(s["exact_allele_accuracy"] - 0.9871) < 1e-12,
            f"synthetic exact accuracy drift: {s['exact_allele_accuracy']}")
    require(s["n_confident"] == 1595, f"synthetic confident-call denominator drift: {s['n_confident']}")
    require(abs(s["false_confident_call_rate"] - 0.0094) < 1e-12,
            f"synthetic false-confident rate drift: {s['false_confident_call_rate']}")

    hifi = read_tsv("results/validation/hifi_short_vs_hifi.tsv")
    compared = [r for r in hifi if r["comparison_universe"] == "YES"
                and r["hifi_truth_status"] == "HIGH_CONFIDENCE_HIFI_TRUTH"]
    exact = [r for r in compared if r["outcome"] == "EXACT_CONCORDANT"]
    require(len(compared) == 913, f"HiFi comparison denominator drift: {len(compared)}")
    require(len(exact) == 911, f"HiFi exact count drift: {len(exact)}")

    catalogue = read_tsv("metadata/loci/repeat_catalogue.tsv")
    require(len(catalogue) == 4009, f"catalogue member count drift: {len(catalogue)}")
    require(len({r["locus_id"] for r in catalogue}) == 1405,
            "catalogue locus count drift")

    audit = read_tsv("metadata/validation/catalogue_boundary_audit.tsv")
    audit_counts = Counter(r["audit_status"] for r in audit)
    expected_audit = {
        "BOUNDARY_STABLE": 400,
        "BOUNDARY_AMBIGUOUS": 48,
        "ORTHOLOGY_AMBIGUOUS": 497,
        "INSUFFICIENT_TO_RESOLVE": 460,
    }
    require(dict(audit_counts) == expected_audit,
            f"boundary-audit counts drift: {dict(audit_counts)}")

    universe = read_tsv("metadata/loci/analysis_universe.tsv")
    universe_counts = Counter(r["g7_locus_class"] for r in universe)
    expected_universe = {
        "PRIMARY_TECHNICAL": 261,
        "PROVISIONAL_TECHNICAL": 39,
        "PRIMARY_INELIGIBLE": 1105,
    }
    require(dict(universe_counts) == expected_universe,
            f"analysis-universe counts drift: {dict(universe_counts)}")
    require(sum(r["known_pv"] == "YES" and r["g7_locus_class"] == "PRIMARY_TECHNICAL"
                for r in universe) == 0,
            "KNOWN_PV loci unexpectedly entered PRIMARY_TECHNICAL")

    screen = read_tsv("results/pilot_disease/locus_screen.tsv")
    passed = Counter((r["cohort"], r["locus_class"]) for r in screen
                     if r["pass_screen"].lower() == "true")
    require(passed[("PRJNA360417", "PRIMARY_TECHNICAL")] == 4,
            f"PRJNA360417 primary screen drift: {passed}")
    require(passed[("PRJNA678459", "PRIMARY_TECHNICAL")] == 6,
            f"PRJNA678459 primary screen drift: {passed}")
    require(passed[("PRJNA1103397", "PRIMARY_TECHNICAL")] == 0,
            f"PRJNA1103397 primary screen drift: {passed}")
    require(passed[("PRJNA360417", "PROVISIONAL_TECHNICAL")] == 0,
            f"PRJNA360417 provisional screen drift: {passed}")
    require(passed[("PRJNA678459", "PROVISIONAL_TECHNICAL")] == 1,
            f"PRJNA678459 provisional screen drift: {passed}")

    cohort = {r["cohort"]: r for r in read_tsv("results/pilot_disease/cohort_level_summary.tsv")}
    require(cohort["PRJNA360417"]["strict_computable"] == "False",
            "PRJNA360417 strict endpoint unexpectedly computable")
    require(cohort["PRJNA678459"]["strict_computable"] == "False",
            "PRJNA678459 strict endpoint unexpectedly computable")

    structure = {r["cohort"]: r for r in read_tsv("results/pilot_disease/structure_overlay_summary.tsv")}
    require(structure["PRJNA678459"]["exact_p"] == "0.031746",
            f"PRJNA678459 structure diagnostic drift: {structure['PRJNA678459']['exact_p']}")
    require(structure["PRJNA360417"]["method"] == "NOT_PERFORMED",
            "PRJNA360417 structure overlay unexpectedly changed")

    print("Publication audit PASS")
    print("synthetic: 1757/1780 exact (98.7%); 15/1595 false-confident (0.94%)")
    print("HiFi: 911/913 exact high-confidence comparisons")
    print("catalogue: 1405 loci / 4009 members")
    print("technical universe: 261 primary / 39 provisional / 1105 ineligible")
    print("pilot primary screens: PRJNA360417=4, PRJNA678459=6, PRJNA1103397=0")


if __name__ == "__main__":
    main()
