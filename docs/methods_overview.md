# Methods overview

## Repeat catalogue

Candidate repeat loci were identified across a reference panel of *H. pylori* genomes and reconciled across orthologous regions. The tracked catalogue, reference flanks, known phase-variable loci, and locus-level eligibility table are stored under `metadata/loci/`.

## Short-read caller

The caller uses sequence anchors flanking each repeat tract to identify spanning fragments and infer tract length. The frozen caller parameters are stored in `config/caller_rules.tsv`. The implementation is in `src/hp_phasevar/`.

Calls require sufficient spanning-read support, acceptable anchor quality, a concentrated allele distribution, and tract lengths compatible with the repeat motif. A noise model derived from within-host data is used to distinguish minor alleles from expected stutter/error.

## Technical validation

Synthetic benchmarking and high-accuracy long-read comparisons were used to assess exact repeat-length concordance. Locus-level technical eligibility additionally incorporates catalogue-boundary stability, orthology resolution, and platform-dependent callability.

The main validation summaries are stored in `results/validation/` and `metadata/validation/`.

## Confounding audit

Disease-labelled public datasets were audited for patient independence, geography, sequencing platform, biopsy site, study membership, and population structure. Summary confounding tables are stored in `results/confounding/`.

## Pilot disease analysis

Analysis rules were frozen before disease-associated testing and are stored in `config/pilot_analysis_rules.tsv`. Only within-study contrasts were considered inferentially interpretable. The tracked pilot outputs are stored in `results/pilot_disease/`.

The disease analysis is exploratory and is not interpreted as evidence for causation or as a definitive pan-population association study.
