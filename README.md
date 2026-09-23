# H. pylori phase-variable repeat analysis

This repository contains the publication-oriented analysis for a study of phase-variable repeat loci in *Helicobacter pylori*.

## Scope

The project develops and validates a short-read framework for calling repeat-length variation at candidate phase-variable loci, evaluates technical callability and platform effects, and applies the framework in a pilot analysis of gastric disease cohorts.

The disease-association component is explicitly exploratory. The available public cohorts were not sufficient for a definitive discovery-and-replication test across the Correa cascade, so the repository does not treat the pilot results as evidence of a causal or general disease association.

## Repository structure

- `src/` — repeat-length caller
- `config/` — frozen caller parameters and analysis rules
- `metadata/` — locus catalogue, eligibility tables, and cohort metadata
- `analysis/` — reproducible analysis scripts
- `results/` — tracked summary results used for the manuscript
- `tests/` — unit tests for the caller
- `docs/` — methods and study notes

Raw sequencing data are not tracked in Git. Public accession identifiers and derived manifests are used to reproduce downloads and analyses.

## Status

The technical framework and pilot analyses are complete. This repository is being prepared as the clean, publication-facing version of the project.
