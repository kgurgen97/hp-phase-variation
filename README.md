# H. pylori repeat-length profiling

This repository contains the publication-oriented code, metadata, and tracked results for a conservative short-read framework to measure repeat-length variation in *Helicobacter pylori*.

## Scope

The project builds a multi-reference repeat catalogue, evaluates read-level callability and study/platform-associated technical variation, validates repeat-length calls against synthetic and high-accuracy long-read evidence, and applies the frozen technical framework in a limited pilot analysis of gastric-disease cohorts.

The catalogue is a **repeat-locus catalogue**, not a claim that every locus is experimentally established as phase variable. The primary technical analysis universe was selected for measurement reliability and should not be interpreted as a set of known phase-variable genes.

The gastric-disease application is exploratory. Available public cohorts did not support a definitive discovery-and-replication test across the Correa cascade, and the repository does not treat the pilot results as evidence of a causal or general disease association.

## Repository structure

- `src/` — repeat-length caller
- `config/` — frozen caller, validation, and pilot-analysis rules
- `metadata/` — repeat catalogue, technical eligibility tables, cohort metadata, and public-data manifests
- `analysis/` — frozen generating-source export, path/provenance maps, and publication audits
- `results/` — tracked manuscript-facing results and the manuscript-number source of truth
- `environment/` — tested publication reproduction baseline
- `tests/` — caller unit tests
- `docs/` — methods and study-scope notes

Raw sequencing data are not tracked. The full Phase-1 call matrix is not duplicated in this clean repository; its exact development-archive blob identity and regeneration route are recorded in `results/calls/phase1_calls.MANIFEST.tsv`. No zero-byte placeholder is retained.

## Quick checks

```bash
python -m unittest discover -s tests
python analysis/check_publication_results.py
python analysis/check_manuscript_numbers.py
```

## Numerical source of truth

Use `results/manuscript/manuscript_numbers.tsv` for the abstract, Results text, tables, and figures. The table is machine-checked against tracked publication artifacts. The 909-locus Phase-1 callable-union entry is explicitly manifest-backed unless the large Phase-1 matrix is materialized.

## Scientific status

The technical framework and limited pilot analysis are frozen for manuscript preparation. The study does **not** establish an association between repeat-length states and gastric carcinogenesis. The full discovery stage was not identifiable from available public cohorts under the pre-specified patient-independence and confounding criteria.
