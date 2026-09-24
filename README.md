# H. pylori repeat-length profiling

This repository contains the publication-oriented code, metadata, and tracked results for a conservative
short-read framework to measure repeat-length variation in *Helicobacter pylori*.

## Scope

The project builds a multi-reference repeat catalogue, evaluates read-level callability and
study/platform-associated technical variation, validates repeat-length calls against synthetic and
high-accuracy long-read evidence, and applies the frozen technical framework in a limited pilot analysis
of gastric-disease cohorts.

The catalogue is a **repeat-locus catalogue**, not a claim that every locus is experimentally established
as phase variable. Likewise, the primary technical analysis universe was selected for measurement
reliability and should not be interpreted as a set of known phase-variable genes.

The gastric-disease application is exploratory. Available public cohorts did not support a definitive
discovery-and-replication test across the Correa cascade, and the repository does not treat the pilot
results as evidence of a causal or general disease association.

## Repository structure

- `src/` — repeat-length caller
- `config/` — frozen caller parameters and pilot-analysis rules
- `metadata/` — repeat catalogue, technical eligibility tables, cohort metadata, and public-data manifests
- `analysis/` — publication consistency/reproducibility scripts
- `results/` — tracked manuscript-facing summary results
- `tests/` — caller unit tests
- `docs/` — methods and study-scope notes

Raw sequencing data and large regenerable call matrices are not tracked in Git. Public accessions,
checksums where available, frozen rules, and the caller implementation are retained. See
`analysis/README.md` for the reproducibility boundary.

## Quick checks

```bash
python -m unittest discover -s tests
python analysis/check_publication_results.py
```

## Scientific status

The technical framework and the limited pilot analysis are frozen for manuscript preparation. The study
does **not** establish an association between repeat-length states and gastric carcinogenesis. The full
discovery stage was not identifiable from the available public cohorts under the pre-specified
patient-independence and confounding criteria.
