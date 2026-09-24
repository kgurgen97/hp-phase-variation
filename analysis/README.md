# Analysis and reproducibility

The publication repository separates frozen generating source, normalized publication-facing inputs/outputs, large regenerable intermediates, and manuscript-number checks.

## Fast checks

```bash
python -m unittest discover -s tests
python analysis/check_publication_results.py
python analysis/check_manuscript_numbers.py
```

These checks do not download reads, change thresholds, redefine disease groups, or compute new biological endpoints.

## Frozen generating source

`analysis/frozen_source/` is a curated export of exact tracked scripts from the authoritative development archive at Project A closeout. It covers catalogue construction, synthetic benchmarking, Phase-1 callability, boundary/orthology audit, HiFi validation, technical locus-universe construction, confounding audit, the bounded disease-blind G7 PRJNA678459 follow-up that generated caller/structure inputs, pilot disease analysis, and structure-diagnostic summarization.

Scripts and the follow-up manifests needed to interpret bounded execution are preserved verbatim. Their internal G5-G8 paths are historical provenance, not a second set of publication results. `analysis/path_map.tsv` maps historical paths to the normalized publication layout, and `analysis/reproduction_manifest.tsv` states the reproduction boundary for each result family.

## Public data manifests

- `metadata/cohorts/phase1_download_manifest.tsv` — Phase-1 short reads
- `metadata/validation/stageA_manifest.tsv` — Stage-A short/ONT validation data
- `metadata/validation/hifi_manifest.tsv` — HiFi/MGI validation data

Raw sequencing files are not versioned.

## Large Phase-1 call matrix

The authoritative historical Phase-1 matrix is `results/g5/g5_phase1_calls.tsv` in the development archive (blob SHA `58486570113e1cde6d0668db1b8d2a6fe38b1915`, 17,870,238 bytes). A zero-byte placeholder is not retained. The exact blob identity and regeneration route are recorded in `results/calls/phase1_calls.MANIFEST.tsv`.

The repository connector could not stream this 17.87 MB blob cross-repository during publication repair. This is a packaging limitation, not a biological result. Manuscript-facing summaries derived from the frozen matrix are tracked.

## Manuscript source of truth

`results/manuscript/manuscript_numbers.tsv` is the machine-readable source for manuscript text, tables, and figures. `analysis/check_manuscript_numbers.py` checks those values against tracked artifacts. The Phase-1 callable-union count is marked as manifest-backed because direct recount requires the omitted large matrix.

## Environment

The exact historical execution environment was only partially recorded and is not reconstructed retrospectively. The newly tested publication baseline is documented under `environment/` and checked in CI on Python 3.10-3.12.

## Interpretation boundary

The catalogue contains repeat loci with different evidence levels. A repeat locus is not automatically an experimentally established phase-variable locus. The `PRIMARY_TECHNICAL` set is a technical analysis universe selected for repeat-length measurement reliability; it contains no literature-anchored `KNOWN_PV` loci. Real-sample outputs should be described as repeat-length alleles, not sample-specific ON/OFF functional states.
