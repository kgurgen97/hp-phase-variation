# Study scope

This project evaluates whether repeat-length variation at technically resolvable repeat loci in *Helicobacter pylori* can be measured reproducibly from public short-read whole-genome sequencing data and explored in gastric-disease cohorts.

The catalogue includes repeat loci with different evidence levels. Repeat-locus membership does not establish phase variation, and the primary technical universe is not an experimentally established phase-variable-gene set.

The publication-facing analysis has two components.

## Technical component

1. Build a multi-reference catalogue of repeat loci.
2. Define read-level anchoring and repeat-length calling rules.
3. Quantify callability and platform/study-associated technical variation.
4. Validate short-read repeat-length calls against high-accuracy long-read sequence where suitable truth data are available.
5. Restrict downstream analyses to loci satisfying the frozen technical eligibility criteria.

## Pilot disease component

The disease analysis is limited to within-cohort contrasts for which the available metadata and sample structure were considered sufficiently interpretable for a pilot stress test. Cross-study pooling was not used for disease inference because disease state was strongly entangled with study, geography, sequencing platform, biopsy site, and/or bacterial population structure.

The pilot analysis is not a definitive test of phase-variable adaptation across the Correa cascade. Under the pre-specified >=10 jointly callable loci per sample-pair floor, the primary cohort-level endpoint was non-computable for both tested contrasts. Publicly available cohorts also did not provide an adequate independent discovery-and-replication design.

Null, non-computable, and inconclusive results are retained as such and were not used to change locus thresholds, disease definitions, or sample inclusion rules.
