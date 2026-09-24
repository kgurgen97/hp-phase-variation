# Main-figure generation

The four main figures are generated from tracked publication artifacts. No script downloads raw reads or changes any frozen biological rule.

Run from the repository root:

```bash
python analysis/figures/figure1_catalogue_universe.py
python analysis/figures/figure2_validation.py
python analysis/figures/figure3_failure_modes.py
python analysis/figures/figure4_public_cohort_stress_test.py
```

Outputs are written as PNG (400 dpi), SVG and PDF under `figures/main/`. Machine-readable panel source tables are written under `figures/source_data/`.

The visual system is defined in `config/figure_style.json`; its BostonGene-associated reference set and limitations are documented in `docs/figure_style_references.md`.

The disease figure is intentionally a stress-test/confounding figure rather than a cancer-association hero figure. Diagnostic p-values are visually secondary and strict non-computable endpoints are not replaced by below-floor calculations.
