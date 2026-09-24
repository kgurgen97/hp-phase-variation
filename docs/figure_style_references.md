# Figure style reference set

The publication figure system is anchored to three BostonGene-associated papers:

1. **Bagaev A et al.** Conserved pan-cancer microenvironment subtypes predict response to immunotherapy. *Cancer Cell*. 2021;39:845–865.e7. DOI: 10.1016/j.ccell.2021.04.014.
2. **Yudina A et al.** Clinical and analytical validation of a combined RNA and DNA exome assay across a large tumor cohort. *Communications Medicine*. 2025;5:236. DOI: 10.1038/s43856-025-00934-3.
3. **Stupichev D et al.** AI-driven multimodal algorithm predicts immunotherapy and targeted therapy outcomes in clear cell renal cell carcinoma. *Cell Reports Medicine*. 2025;6:102299. DOI: 10.1016/j.xcrm.2025.102299.

The useful visual features are consistent across this reference set: white backgrounds, thin axes, disciplined typography, compact multi-panel assembly, muted categorical colors, direct denominators/annotations, and saturated color reserved for meaningful contrasts.

The palette in `config/figure_style.json` is **reference-derived from published raster figures** (especially Yudina Figs. 1–2 and the Bagaev/Stupichev heatmap/annotation panels). Because the published figures are rasterized/anti-aliased and no single official HEX dictionary is provided for all panels, the palette is not described as an official BostonGene brand palette.

## Project-specific constraints

- PRIMARY_TECHNICAL is blue; PROVISIONAL_TECHNICAL is teal; ineligible/null is gray.
- Evidence class is visually separate from technical eligibility: REPEAT_ONLY gray, STRONG_PV_CANDIDATE magenta, KNOWN_PV gold.
- Boundary/orthology failure classes use coral/purple/gray, while BOUNDARY_STABLE uses blue.
- HiFi reference-length comparisons are light blue; non-reference comparisons are magenta.
- Diagnostic disease-association results remain visually subdued and must not dominate the figure hierarchy.
- Population-structure association in PRJNA678459 is a confounding caution flag, not a disease finding.
