#!/usr/bin/env python3
"""Validate manuscript_numbers.tsv against frozen tracked publication artifacts.

Consistency audit only: no downloads, threshold changes, disease reclassification,
or new biological analyses.
"""
from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def read_tsv(path: str):
    with (ROOT / path).open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    numbers = {r["metric"]: r for r in read_tsv("results/manuscript/manuscript_numbers.tsv")}
    failures = []

    def expect(metric: str, value):
        got = numbers[metric]["value"]
        if isinstance(value, bool):
            want = "TRUE" if value else "FALSE"
        elif isinstance(value, float):
            want = format(value, ".12g")
        else:
            want = str(value)
        if got != want:
            failures.append(f"{metric}: manuscript={got!r}, computed={want!r}")

    panel = read_tsv("metadata/loci/reference_panel.tsv")
    expect("reference_panel_genomes", len(panel))
    if len({r["short_code"] for r in panel}) != len(panel):
        failures.append("reference panel contains duplicate short_code values")

    catalogue = read_tsv("metadata/loci/repeat_catalogue.tsv")
    by_locus = defaultdict(list)
    for row in catalogue:
        by_locus[row["locus_id"]].append(row)
    expect("catalogue_members", len(catalogue))
    expect("catalogue_loci", len(by_locus))
    expect("caller_eligible_loci", sum(
        any((r.get("caller_eligible") or "").startswith("YES") for r in members)
        for members in by_locus.values()
    ))

    universe = read_tsv("metadata/loci/analysis_universe.tsv")
    classes = Counter(r["g7_locus_class"] for r in universe)
    expect("primary_technical_loci", classes["PRIMARY_TECHNICAL"])
    expect("provisional_technical_loci", classes["PROVISIONAL_TECHNICAL"])
    expect("primary_ineligible_loci", classes["PRIMARY_INELIGIBLE"])
    expect("primary_known_pv_loci", sum(
        r["g7_locus_class"] == "PRIMARY_TECHNICAL" and r["known_pv"] == "YES"
        for r in universe
    ))
    primary = {r["locus_id"] for r in universe if r["g7_locus_class"] == "PRIMARY_TECHNICAL"}
    evidence = Counter()
    for locus_id in primary:
        member_classes = {r["evidence_class"] for r in by_locus[locus_id]}
        evidence["MIXED" if len(member_classes) > 1 else next(iter(member_classes))] += 1
    expect("primary_repeat_only_loci", evidence["REPEAT_ONLY"])
    expect("primary_strong_pv_candidate_loci", evidence["STRONG_PV_CANDIDATE"])
    expect("primary_mixed_evidence_loci", evidence["MIXED"])

    audit = Counter(r["audit_status"] for r in read_tsv("metadata/validation/catalogue_boundary_audit.tsv"))
    expect("boundary_stable_loci", audit["BOUNDARY_STABLE"])
    expect("boundary_ambiguous_loci", audit["BOUNDARY_AMBIGUOUS"])
    expect("orthology_ambiguous_loci", audit["ORTHOLOGY_AMBIGUOUS"])
    expect("insufficient_to_resolve_loci", audit["INSUFFICIENT_TO_RESOLVE"])

    synthetic = {r["metric"]: r for r in read_tsv("results/validation/synthetic_authoritative_counts.tsv")}
    expect("synthetic_empirical_total", synthetic["empirical_total"]["value"])
    expect("synthetic_callable", synthetic["callable"]["value"])
    expect("synthetic_exact_dominant", synthetic["exact_dominant"]["numerator"])
    expect("synthetic_exact_dominant_fraction", synthetic["exact_dominant"]["value"])
    expect("synthetic_confident", synthetic["confident"]["value"])
    expect("synthetic_false_confident_wrong", synthetic["false_confident_wrong"]["numerator"])
    expect("synthetic_false_confident_wrong_fraction", synthetic["false_confident_wrong"]["value"])
    expect("synthetic_uncallable_rate", synthetic["uncallable_rate"]["value"])
    sj = json.loads((ROOT / "results/validation/synthetic_overall.json").read_text())
    expect("synthetic_uncallable_rate", sj["all"]["uncallable_rate"])

    hifi_j = json.loads((ROOT / "results/validation/hifi_overall.json").read_text())
    truth = hifi_j["truth_status_counts"]
    expect("hifi_truth_high_confidence", truth["HIGH_CONFIDENCE_HIFI_TRUTH"])
    expect("hifi_truth_ambiguous", truth["AMBIGUOUS_HIFI"])
    expect("hifi_truth_uncallable", truth["UNCALLABLE_HIFI"])
    expect("hifi_truth_total", sum(truth.values()))

    hifi = read_tsv("results/validation/hifi_short_vs_hifi.tsv")
    compared = [r for r in hifi
                if r["comparison_universe"] == "YES"
                and r["hifi_truth_status"] == "HIGH_CONFIDENCE_HIFI_TRUTH"]
    expect("hifi_high_confidence_comparisons", len(compared))
    expect("hifi_exact", sum(r["outcome"] == "EXACT_CONCORDANT" for r in compared))
    nr = {r["stratum"]: r for r in read_tsv("results/validation/hifi_nonreference_summary.tsv")}
    nr_all = nr["ALL_HIGH_CONFIDENCE (not a validation of any single class)"]
    expect("hifi_nonreference_comparisons", int(nr_all["nonreference_comparisons"]))
    expect("hifi_nonreference_exact", int(nr_all["nonreference_exact"]))

    callability = {r["metric"]: r for r in read_tsv("results/callability/summary.tsv")}
    expect("phase1_callable_union_loci", callability["phase1_callable_union_loci"]["value"])
    platform = read_tsv("results/callability/platform_transfer.tsv")
    expect("platform_study_confounded_callability_loci",
           sum(r["callability_class"] == "PLATFORM_CONFOUNDED_CALLABILITY" for r in platform))

    conf = read_tsv("results/confounding/confounding_matrix.tsv")
    expect("confounding_correa_studies", len(conf))
    expect("confounding_correa_biosamples", sum(int(r["n_correa_biosamples"]) for r in conf))
    single_country = sum(";" not in r["countries"] for r in conf)
    expect("confounding_single_country_studies", single_country)
    expect("confounding_multicountry_studies", len(conf) - single_country)

    pilot = json.loads((ROOT / "results/pilot_disease/summary.json").read_text())
    expect("prjna360417_nag", pilot["prjna360417"]["n_nag"])
    expect("prjna360417_im", pilot["prjna360417"]["n_im"])
    expect("prjna678459_ag", pilot["prjna678459"]["n_ag"])
    expect("prjna678459_gc", pilot["prjna678459"]["n_gc"])
    expect("prjna1103397_gc", pilot["prjna1103397"]["n_gc"])

    cohort = {r["cohort"]: r for r in read_tsv("results/pilot_disease/cohort_level_summary.tsv")}
    c360, c678, c110 = cohort["PRJNA360417"], cohort["PRJNA678459"], cohort["PRJNA1103397"]
    expect("prjna360417_primary_screen_pass", int(c360["n_primary_passed_screen"]))
    expect("prjna360417_strict_computable", c360["strict_computable"].upper())
    expect("prjna360417_diagnostic_exact_p", float(c360["diagnostic_exact_p"]))
    expect("prjna678459_primary_screen_pass", int(c678["n_primary_passed_screen"]))
    expect("prjna678459_provisional_screen_pass", int(c678["n_provisional_passed_screen"]))
    expect("prjna678459_strict_computable", c678["strict_computable"].upper())
    expect("prjna678459_diagnostic_exact_p", float(c678["diagnostic_exact_p"]))
    expect("prjna1103397_primary_screen_pass", int(c110["n_primary_passed_screen"]))

    structure = {r["cohort"]: r for r in read_tsv("results/pilot_disease/structure_overlay_summary.tsv")}
    s = structure["PRJNA678459"]
    expect("structure_prjna678459_pseudo_f", float(s["pseudo_F"]))
    expect("structure_prjna678459_r2", float(s["R2"]))
    expect("structure_prjna678459_exact_p", float(s["exact_p"]))
    expect("structure_prjna678459_permutations", int(s["n_permutations"]))
    expect("structure_prjna360417", structure["PRJNA360417"]["method"])

    if failures:
        raise AssertionError("\n".join(failures))

    print(f"Manuscript number audit PASS ({len(numbers)} frozen entries)")
    print("Phase-1 callable-union count is manifest-backed unless the 17.87 MB call matrix is materialized.")


if __name__ == "__main__":
    main()
