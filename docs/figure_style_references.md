# Figure style reference set

The publication figure system is anchored to three BostonGene-associated papers:

1. **Bagaev A et al.** *Conserved pan-cancer microenvironment subtypes predict response to immunotherapy*. Cancer Cell. 2021;39:845–865.e7. DOI: 10.1016/j.ccell.2021.04.014.
2. **Yudina A et al.** *Clinical and analytical validation of a combined RNA and DNA exome assay across a large tumor cohort*. Communications Medicine. 2025;5:236. DOI: 10.1038/s43856-025-00934-3.
3. **Stupichev D et al.** *AI-driven multimodal algorithm predicts immunotherapy and targeted therapy outcomes in clear cell renal cell carcinoma*. Cell Reports Medicine. 2025;6:102299. DOI: 10.1016/j.xcrm.2025.102299.

## Exact palette provenance

The project palette is **not estimated from raster screenshots**. Exact HEX values in `config/figure_style.json` were parsed as literal `#RRGGBB` tokens from BostonGene's public MFP companion SVG, `BostonGene/MFP:img/Abstract.svg`, blob SHA `b01665d37017c0e69f3703ef096845925e93ba38`, associated with Bagaev et al. The dominant non-black/non-white literals include `#333F50`, `#5555CC`, `#36BDC0`, `#CA952A`, `#CA2C50`, `#EBA382`, `#2F6354`, `#DE6A36`, `#863541`, `#56D3BB`, `#71C293`, `#D689B1`, and `#84ACE8`.

Yudina 2025 and Stupichev 2025 are used as layout references rather than as sources of guessed HEX values. Their published figures support the same design grammar: white backgrounds, thin axes, compact multi-panel assembly, direct denominator annotations, disciplined heatmaps, and saturated accents reserved for meaningful contrasts.

## Project-specific constraints

- PRIMARY_TECHNICAL uses indigo; PROVISIONAL_TECHNICAL uses teal; ineligible/null uses light gray.
- Evidence class remains visually separate from technical eligibility: REPEAT_ONLY gray, STRONG_PV_CANDIDATE mauve, KNOWN_PV gold.
- Boundary/orthology failure classes use orange/plum/gray; resolved technical states use teal/green.
- HiFi reference-length and non-reference comparisons use distinct colors and always show denominators.
- Diagnostic disease-association results remain visually subdued and must not dominate the figure hierarchy.
- The PRJNA678459 structure statistic is a confounding caution flag, not a repeat-disease finding.
- PRJNA360417 structure remains NOT_PERFORMED.

The plotting code performs descriptive aggregation for visualization only. It does not change cohorts, thresholds, disease definitions, endpoints, or frozen inferential results.
