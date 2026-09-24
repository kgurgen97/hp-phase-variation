# Methods overview

## Repeat catalogue

Repeat loci were identified across a reference panel of *H. pylori* genomes and reconciled across orthologous regions. Evidence annotations distinguish repeat-only loci, stronger phase-variation candidates, and literature-anchored known phase-variable loci. The tracked catalogue, reference flanks, known-PV annotations, and locus-level eligibility tables are stored under `metadata/loci/`.

## Short-read caller

The caller uses sequence anchors flanking each repeat tract to identify spanning fragments and infer repeat length. Frozen caller parameters are stored in `config/caller_rules.tsv`, and the implementation is in `src/hp_phasevar/`.

Calls require sufficient spanning-read support, acceptable anchor quality, a concentrated allele distribution, and tract lengths compatible with the repeat motif. A noise model derived during development is used to distinguish minor alleles from expected stutter/error. Real-sample outputs are interpreted as repeat-length alleles; sample-specific ON/OFF functional states are not claimed.

## Technical validation

Synthetic benchmarking and high-accuracy long-read comparisons assess exact repeat-length concordance. Locus-level technical eligibility additionally incorporates catalogue-boundary stability, orthology resolution, and platform/study-confounded callability. The latter is descriptive technical association and is not interpreted as a causal platform effect.

Validation summaries are stored in `results/validation/`, `results/callability/`, and `metadata/validation/`. Exact frozen generating-source exports are under `analysis/frozen_source/`.

## Confounding audit

Disease-labelled public datasets were audited for patient independence, geography, sequencing platform, biopsy site, study membership, and screening-level population structure. The MinHash/k-mer structure analysis is a confounding screen, not an ancestry-adjusted disease model. Cross-study disease inference was not performed.

## Pilot disease analysis

Analysis rules were frozen before disease-associated testing and are stored in `config/pilot_analysis_rules.tsv`. The measurable object is the repeat-length allele, and the cohort-level distance is described as repeat-allele-state or repeat-length allele distance.

The pre-specified >=10 jointly callable loci per sample-pair floor gated the strict cohort-level endpoint. Both pilot contrasts failed that floor, so their strict endpoints are reported as non-computable; below-floor PERMANOVA values are retained only as diagnostics. PRJNA360417 has no population-structure overlay and remains an explicit limitation.

The pilot analysis is exploratory and is not interpreted as evidence for causation or as a definitive pan-population association study.
