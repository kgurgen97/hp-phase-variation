# Publication reproduction environment

## Historical execution environment

The development-stage environment was only partially recorded. Exact historical Python patch versions and a complete historical library lockfile are therefore not claimed retrospectively.

## Newly tested publication baseline

The publication checks are tested in GitHub Actions on Python 3.10, 3.11, and 3.12. Python 3.11 is the declared publication baseline (see `.python-version` and `environment/python-version.txt`).

The caller unit tests and publication consistency audits currently use only the Python standard library. This baseline is a newly tested publication reproduction environment, not a reconstruction of every historical development-gate environment.

From-raw regeneration is stage-specific and may require external command-line tools plus public sequencing/reference data. Accessions, checksums where available, frozen rules, and the exact generating source are tracked. Publication CI performs no raw-read download and no new biological analysis.
