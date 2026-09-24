# Analysis and reproducibility

The publication repository separates **tracked manuscript-facing outputs** from large, regenerable
intermediate call matrices.

## Fast consistency check

Run:

```bash
python analysis/check_publication_results.py
python -m unittest discover -s tests
```

The audit verifies the key numerical claims used in the manuscript against the tracked TSV/JSON files.
It does not download reads or recompute disease associations.

## Public data manifests

The repository tracks the accessions and download metadata needed for the major raw-data stages:

- `metadata/cohorts/phase1_download_manifest.tsv` — Phase-1 short reads;
- `metadata/validation/stageA_manifest.tsv` — Stage-A short/ONT validation data;
- `metadata/validation/hifi_manifest.tsv` — HiFi/MGI validation data.

Raw sequencing files are intentionally not versioned.

## Large derived matrices

The full Phase-1 per-locus call matrix is a large, regenerable intermediate and is not stored in Git.
It is produced by the frozen caller from the public reads and frozen catalogue/rules. Manuscript-facing
summary outputs derived from the frozen analysis are tracked under `results/`.

The repository does not claim that the current environment reconstructs undocumented historical software
versions used during every development gate. The code and tracked results define the publication release;
historical environment gaps are disclosed in the manuscript/reproducibility notes.

## Interpretation boundary

The catalogue contains repeat loci with different evidence levels. A repeat locus is not automatically an
experimentally established phase-variable locus. The `PRIMARY_TECHNICAL` set is a technical analysis
universe selected for repeat-length measurement reliability; it contains no literature-anchored
`KNOWN_PV` loci.
