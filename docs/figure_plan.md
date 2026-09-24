# Frozen publication figure plan

The figure set is designed for the methods/resource framing of Project A. It does not present the pilot disease analysis as a positive cancer-association result.

## Figure 1 — Catalogue, audit, and technical universe

**A. Catalogue scale.** 1,405 repeat loci / 4,009 catalogue members; 1,196 loci with at least one caller-eligible member.

**B. Boundary/orthology audit.** Four frozen audit states across the full catalogue.

**C. Frozen technical universe.** PRIMARY_TECHNICAL, PROVISIONAL_TECHNICAL, and PRIMARY_INELIGIBLE counts.

**D. Evidence composition within PRIMARY_TECHNICAL.** REPEAT_ONLY, STRONG_PV_CANDIDATE, mixed member-level evidence, and the explicit absence of literature-anchored KNOWN_PV loci from PRIMARY.

The layout deliberately separates evidence level from technical eligibility so the figure cannot imply that PRIMARY_TECHNICAL is an experimentally established phase-variable gene set.

## Figure 2 — Technical validation

**A. Synthetic empirical benchmark.** Callable fraction, exact dominant repeat-length concordance among callable comparisons, confidence fraction among callable comparisons, and false-confident error among confident calls. Numerators and denominators are printed on the panel.

**B. HiFi truth availability.** HIGH_CONFIDENCE_HIFI_TRUTH, AMBIGUOUS_HIFI, and UNCALLABLE_HIFI across all isolate-locus pairs.

**C. Short-read versus HiFi concordance.** Overall exact concordance and non-reference exact concordance, with denominators.

**D. Non-reference validation coverage by motif stratum.** This panel distinguishes strata with observed non-reference truth from strata for which non-reference alleles were not validated.

## Figure 3 — Technical failure modes and interpretation boundary

**A. Boundary/orthology failure-state distribution.**

**B. Study/platform-associated callability classes.** Counts are taken from the frozen platform-transfer table and are descriptive technical associations, not causal platform effects.

**C. Known-PV representation across frozen technical classes.** This makes the evidence gap visible rather than obscuring it.

**D. Reproducibility boundary.** A concise text panel states that the exact historical environment is incomplete, the large Phase-1 matrix is manifest-backed, and G6 curated closure-table provenance remains disclosed.

## Figure 4 — Public-cohort stress test

**A. Confounding landscape.** 45 Correa-labelled studies, 793 BioSamples, with single-country versus multicountry study counts.

**B. Pilot cohort sizes and available contrasts.**

**C. Locus-screen attrition.** PRIMARY screen-pass counts for each pilot cohort, plus the one PROVISIONAL locus in PRJNA678459.

**D. Endpoint status.** Both tested within-cohort contrasts remain NOT_COMPUTABLE under the pre-specified >=10 jointly callable loci per sample-pair floor. Below-floor PERMANOVA values are displayed only as diagnostics; PRJNA678459 structure is a confounding caution flag; PRJNA360417 structure remains NOT_PERFORMED.

## Visual reference provenance

The style configuration uses exact categorical colors reproduced from the public BostonGene/MFP plotting implementation at commit `4debc39bb1bae10dfc550ec9fcb1c653d722c773`, where the default categorical palette is generated from `matplotlib.cm.rainbow` over `[0, 0.92]` and clustered heatmaps use `coolwarm`.

Published figure layouts were also reviewed from:

- Bagaev et al. *Cancer Cell* (2021), doi:10.1016/j.ccell.2021.04.014.
- Yudina et al. *Communications Medicine* (2025), doi:10.1038/s43856-025-00934-3.
- Stupichev et al. *Cell Reports Medicine* (2025), doi:10.1016/j.xcrm.2025.102299.

These papers are visual references only; no scientific inference or data are imported from them.
