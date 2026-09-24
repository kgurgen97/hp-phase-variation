#!/usr/bin/env python3
"""G7 section 4: disease-blind genome-wide k-mer/sketch distance analysis.

Used only because no established whole-genome/core-genome tool (mash,
sourmash, skani, fastANI, parsnp, etc.) is installed in this environment and
no packages may be installed in this gate (G7 issue, preferred-hierarchy
option C). Labelled SCREENING_STRUCTURE_ONLY throughout; this is NOT a
core-genome SNP phylogeny and must not be read as clonal ancestry (H. pylori
is highly recombinogenic).

Method (own implementation, no external bioinformatics package, no numpy):
- Each assembly FASTA (gzip) is read; contigs are kept separate so k-mers
  never cross a contig boundary.
- For every k-mer (k=KMER_SIZE) not containing N, the canonical form is the
  lexicographic min of the k-mer and its reverse complement.
- A bottom-SKETCH_SIZE MinHash sketch (the SKETCH_SIZE smallest zlib.crc32
  hash values of the canonical k-mers) is kept per genome (deterministic
  hash, not Python's randomized hash()).
- Pairwise distance = 1 - (Jaccard estimate from sketch intersection/union).

Run: python3 scripts/g7/g7_kmer_sketch_distance.py
"""
import csv
import gzip
import heapq
import json
import os
import time
import zlib

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
KMER_SIZE = 21
SKETCH_SIZE = 1000
RC_TABLE = str.maketrans("ACGT", "TGCA")


def read_fasta_gz_sequences(path):
    seqs = []
    cur = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur))
                    cur = []
            else:
                cur.append(line.strip().upper())
    if cur:
        seqs.append("".join(cur))
    return seqs


def sketch_genome(seqs, k=KMER_SIZE, sketch_size=SKETCH_SIZE):
    seen = set()
    max_heap = []
    for seq in seqs:
        n = len(seq)
        if n < k:
            continue
        for i in range(0, n - k + 1):
            kmer = seq[i:i + k]
            if "N" in kmer:
                continue
            rc = kmer.translate(RC_TABLE)[::-1]
            canon = kmer if kmer < rc else rc
            hv = zlib.crc32(canon.encode("ascii"))
            if hv in seen:
                continue
            if len(max_heap) < sketch_size:
                heapq.heappush(max_heap, -hv)
                seen.add(hv)
            elif -hv > max_heap[0]:
                removed = -heapq.heapreplace(max_heap, -hv)
                seen.discard(removed)
                seen.add(hv)
    return seen


def jaccard(a, b):
    if not a or not b:
        return 0.0
    inter = len(a & b)
    union = len(a | b)
    return inter / union if union else 0.0


def main():
    manifest_path = os.path.join(REPO, "metadata", "g7", "g7_population_structure_manifest.tsv")
    with open(manifest_path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    rows = [r for r in rows if r.get("local_fasta_path") and r["local_fasta_path"] != "NA"]

    sketches = {}
    timing = {}
    for r in rows:
        sid = r["sample_id"]
        path = os.path.join(REPO, r["local_fasta_path"])
        t0 = time.time()
        seqs = read_fasta_gz_sequences(path)
        sk = sketch_genome(seqs, k=KMER_SIZE, sketch_size=SKETCH_SIZE)
        timing[sid] = round(time.time() - t0, 2)
        sketches[sid] = sk
        print(sid, "sketch_size", len(sk), "seconds", timing[sid], flush=True)

    ids = list(sketches.keys())
    out_matrix = os.path.join(REPO, "results", "g7", "g7_distance_matrix.tsv")
    with open(out_matrix, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample_id"] + ids)
        for i in ids:
            row = [i]
            for j in ids:
                d = 1.0 - jaccard(sketches[i], sketches[j])
                row.append(f"{d:.6f}")
            w.writerow(row)

    timing_path = os.path.join(REPO, "results", "g7", "g7_kmer_sketch_timing.json")
    with open(timing_path, "w") as fh:
        json.dump(
            {"kmer_size": KMER_SIZE, "sketch_size": SKETCH_SIZE, "hash": "zlib.crc32",
             "method": "SCREENING_STRUCTURE_ONLY_bottom_minhash_own_implementation",
             "seconds_per_genome": timing},
            fh, indent=2)

    print("done", len(ids), "genomes")


if __name__ == "__main__":
    main()
