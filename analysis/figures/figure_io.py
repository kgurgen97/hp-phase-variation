from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def read_tsv(path):
    with (ROOT / path).open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def manuscript_numbers():
    return {r["metric"]: r for r in read_tsv("results/manuscript/manuscript_numbers.tsv")}


def n_int(numbers, metric):
    return int(numbers[metric]["value"])


def n_float(numbers, metric):
    return float(numbers[metric]["value"])


def write_tsv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def count_primary_motifs():
    rows = read_tsv("metadata/loci/analysis_universe.tsv")
    return Counter(r["motif_class"] for r in rows if r["g7_locus_class"] == "PRIMARY_TECHNICAL")
