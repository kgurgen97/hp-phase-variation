# Publication figures

This directory contains the publication-facing rendering code and its visual style configuration.

The scripts must not recompute biological endpoints. They read frozen tracked outputs and the machine-checked manuscript-number table, then render figures from those values.

## Render

```bash
python -m pip install -r environment/figure-requirements.txt
python analysis/figures/make_publication_figures.py --output-dir figures/generated
```

The script writes SVG, PDF, and PNG versions of Figures 1–4.

## Numerical contract

Primary counts, denominators, and pilot statistics are read from `results/manuscript/manuscript_numbers.tsv`. Supporting categorical distributions are read directly from tracked frozen TSVs. Figure rendering must fail loudly if a required field is absent.

No threshold, disease label, cohort inclusion decision, or statistical endpoint is changed by this code.
