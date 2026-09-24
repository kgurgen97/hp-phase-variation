#!/usr/bin/env python3
from __future__ import annotations

import csv
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/"tables"/"main"
OUT.mkdir(parents=True,exist_ok=True)

def read(path):
    with (ROOT/path).open(newline="") as fh:
        return list(csv.DictReader(fh,delimiter="\t"))

def write(name,rows,fields):
    with (OUT/name).open("w",newline="") as fh:
        w=csv.DictWriter(fh,fields,delimiter="\t",lineterminator="\n")
        w.writeheader();w.writerows(rows)

N={r["metric"]:r for r in read("results/manuscript/manuscript_numbers.tsv")}

def val(k): return N[k]["value"]

table1=[
 {"section":"Catalogue","metric":"Repeat loci","value":val("catalogue_loci"),"denominator":"","note":"Multi-reference repeat-locus catalogue"},
 {"section":"Catalogue","metric":"Catalogue members","value":val("catalogue_members"),"denominator":"","note":"Reference-panel members"},
 {"section":"Catalogue","metric":"Caller-eligible loci","value":val("caller_eligible_loci"),"denominator":"","note":"Initial G5 technical eligibility"},
 {"section":"Boundary/orthology","metric":"BOUNDARY_STABLE","value":val("boundary_stable_loci"),"denominator":val("catalogue_loci"),"note":"Disease-blind audit"},
 {"section":"Boundary/orthology","metric":"BOUNDARY_AMBIGUOUS","value":val("boundary_ambiguous_loci"),"denominator":val("catalogue_loci"),"note":"Disease-blind audit"},
 {"section":"Boundary/orthology","metric":"ORTHOLOGY_AMBIGUOUS","value":val("orthology_ambiguous_loci"),"denominator":val("catalogue_loci"),"note":"Disease-blind audit"},
 {"section":"Boundary/orthology","metric":"INSUFFICIENT_TO_RESOLVE","value":val("insufficient_to_resolve_loci"),"denominator":val("catalogue_loci"),"note":"Disease-blind audit"},
 {"section":"Final universe","metric":"PRIMARY_TECHNICAL","value":val("primary_technical_loci"),"denominator":val("catalogue_loci"),"note":"Primary technical repeat-locus universe"},
 {"section":"Final universe","metric":"PROVISIONAL_TECHNICAL","value":val("provisional_technical_loci"),"denominator":val("catalogue_loci"),"note":"Prespecified sensitivity only"},
 {"section":"Final universe","metric":"PRIMARY_INELIGIBLE","value":val("primary_ineligible_loci"),"denominator":val("catalogue_loci"),"note":"Excluded technically"},
 {"section":"Synthetic validation","metric":"Exact dominant allele","value":val("synthetic_exact_dominant"),"denominator":val("synthetic_callable"),"note":"98.7% of callable empirical comparisons"},
 {"section":"Synthetic validation","metric":"False-confident wrong dominant allele","value":val("synthetic_false_confident_wrong"),"denominator":val("synthetic_confident"),"note":"0.94% of HIGH/MEDIUM-confidence calls"},
 {"section":"HiFi validation","metric":"Exact high-confidence comparisons","value":val("hifi_exact"),"denominator":val("hifi_high_confidence_comparisons"),"note":"Descriptive per-comparison concordance"},
 {"section":"HiFi validation","metric":"Exact non-reference comparisons","value":val("hifi_nonreference_exact"),"denominator":val("hifi_nonreference_comparisons"),"note":"Non-reference repeat lengths"},
]
write("table1_technical_summary.tsv",table1,["section","metric","value","denominator","note"])

cohort={r["cohort"]:r for r in read("results/pilot_disease/cohort_level_summary.tsv")}
structure={r["cohort"]:r for r in read("results/pilot_disease/structure_overlay_summary.tsv")}
table2=[
 {"cohort":"PRJNA360417","contrast":"NAG vs IM","group_sizes":"5 vs 6","primary_passed":cohort["PRJNA360417"]["n_primary_passed_screen"],
  "provisional_passed":cohort["PRJNA360417"]["n_provisional_passed_screen"],"pairwise_coverage":f"{cohort['PRJNA360417']['min_pairwise_coverage']}-{cohort['PRJNA360417']['max_pairwise_coverage']}",
  "strict_endpoint":"NOT_COMPUTABLE","diagnostic_exact_p":cohort["PRJNA360417"]["diagnostic_exact_p"],
  "structure_status":structure["PRJNA360417"]["method"],"interpretation":"Below prespecified >=10-locus pairwise floor"},
 {"cohort":"PRJNA678459","contrast":"AG vs GC","group_sizes":"5 vs 5","primary_passed":cohort["PRJNA678459"]["n_primary_passed_screen"],
  "provisional_passed":cohort["PRJNA678459"]["n_provisional_passed_screen"],"pairwise_coverage":f"{cohort['PRJNA678459']['min_pairwise_coverage']}-{cohort['PRJNA678459']['max_pairwise_coverage']}",
  "strict_endpoint":"NOT_COMPUTABLE","diagnostic_exact_p":cohort["PRJNA678459"]["diagnostic_exact_p"],
  "structure_status":"SCREENING_STRUCTURE_ONLY; exact p="+structure["PRJNA678459"]["exact_p"],
  "interpretation":"Below prespecified >=10-locus pairwise floor; structure result is a confounding caution"},
 {"cohort":"PRJNA1103397","contrast":"GC only","group_sizes":"6 isolates; one selected isolate per patient reported","primary_passed":cohort["PRJNA1103397"]["n_primary_passed_screen"],
  "provisional_passed":"NA","pairwise_coverage":"NA","strict_endpoint":"NO_CONTRAST","diagnostic_exact_p":"NA",
  "structure_status":"NOT_APPLICABLE","interpretation":"Explicit patient identifiers unavailable; no disease contrast"},
]
write("table2_pilot_summary.tsv",table2,["cohort","contrast","group_sizes","primary_passed","provisional_passed","pairwise_coverage","strict_endpoint","diagnostic_exact_p","structure_status","interpretation"])

table3=[
 {"domain":"Public Correa-labelled data","metric":"Studies","value":val("confounding_correa_studies"),"note":"Audited studies"},
 {"domain":"Public Correa-labelled data","metric":"BioSamples","value":val("confounding_correa_biosamples"),"note":"Audited BioSamples"},
 {"domain":"Geography","metric":"Single-country studies","value":f"{val('confounding_single_country_studies')}/{val('confounding_correa_studies')}","note":"Current tracked table; supersedes stale 44/45 handoff value"},
 {"domain":"Geography","metric":"Multi-country studies","value":f"{val('confounding_multicountry_studies')}/{val('confounding_correa_studies')}","note":"Four studies"},
 {"domain":"Phase-1","metric":"Callable union","value":val("phase1_callable_union_loci"),"note":"Manifest-backed unless full Phase-1 matrix materialized"},
 {"domain":"Phase-1","metric":"Platform/study-confounded callability","value":val("platform_study_confounded_callability_loci"),"note":"Descriptive technical association, not causal platform effect"},
 {"domain":"Discovery/replication","metric":"Adequate G9 discovery set","value":"NONE","note":"BLOCKED_NO_ADEQUATE_DISCOVERY_DATA"},
 {"domain":"Discovery/replication","metric":"Independent replication set","value":"NONE","note":"No replication analysis performed"},
]
write("table3_data_limits.tsv",table3,["domain","metric","value","note"])

print("Wrote",*(p.name for p in sorted(OUT.glob("table*.tsv"))))
