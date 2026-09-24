#!/usr/bin/env python3
"""G7 follow-up: disease-blind population-structure update. Reuses the frozen G7 SCREENING_STRUCTURE_ONLY
k-mer/MinHash sketch method (scripts/g7/g7_kmer_sketch_distance.py, unchanged) to place the 10 PRJNA678459
samples (read-based sketches, bounded to 30000 read pairs/sample; see g7f_scope_rules.tsv F07) relative to
the existing 32-genome G7 distance matrix. No disease label is read anywhere in this script (F06)."""
import csv, gzip, heapq, json, os, sys, time, zlib

sys.path.insert(0, 'scripts/g7')
from g7_kmer_sketch_distance import sketch_genome, jaccard, KMER_SIZE, SKETCH_SIZE  # noqa: E402

RAW = 'data/raw/g7_followup/PRJNA678459'
READ_PAIRS_CAP = 30000


def read_fastq_seqs(path, max_records):
    seqs = []
    with gzip.open(path, 'rt') as fh:
        n = 0
        while n < max_records:
            h = fh.readline()
            if not h:
                break
            s = fh.readline().rstrip('\n')
            fh.readline()
            fh.readline()
            seqs.append(s)
            n += 1
    return seqs


def main():
    qc_rows = list(csv.DictReader(open('metadata/g7_followup/g7f_qc_manifest.tsv'), delimiter='\t'))

    # rebuild the 32 existing genome sketches (unchanged method; not persisted by the original G7 run)
    manifest = list(csv.DictReader(open('metadata/g7/g7_population_structure_manifest.tsv'), delimiter='\t'))
    manifest = [r for r in manifest if r.get('local_fasta_path') and r['local_fasta_path'] != 'NA']
    sketches = {}
    for r in manifest:
        sid = r['sample_id']
        from g7_kmer_sketch_distance import read_fasta_gz_sequences
        seqs = read_fasta_gz_sequences(r['local_fasta_path'])
        sketches[sid] = sketch_genome(seqs, k=KMER_SIZE, sketch_size=SKETCH_SIZE)
        print('genome', sid, len(sketches[sid]), flush=True)

    timing = {}
    for r in qc_rows:
        biosample, run = r['biosample'], r['run_accession']
        r1 = os.path.join(RAW, run, run + '_1.fastq.gz')
        r2 = os.path.join(RAW, run, run + '_2.fastq.gz')
        t0 = time.time()
        seqs = read_fastq_seqs(r1, READ_PAIRS_CAP) + read_fastq_seqs(r2, READ_PAIRS_CAP)
        sk = sketch_genome(seqs, k=KMER_SIZE, sketch_size=SKETCH_SIZE)
        timing[biosample] = round(time.time() - t0, 2)
        sketches[biosample] = sk
        print('read_sample', biosample, len(sk), timing[biosample], 's', flush=True)

    genome_ids = [r['sample_id'] for r in manifest]
    new_ids = [r['biosample'] for r in qc_rows]
    all_ids = genome_ids + new_ids

    with open('results/g7_followup/g7f_structure_distance_matrix.tsv', 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['sample_id'] + all_ids)
        for i in all_ids:
            row = [i] + [f'{1.0 - jaccard(sketches[i], sketches[j]):.6f}' for j in all_ids]
            w.writerow(row)

    # nearest-neighbor placement of each new read-based sample among the 32 existing genome sketches only
    nn_rows = []
    for nid in new_ids:
        dists = sorted(((1.0 - jaccard(sketches[nid], sketches[gid]), gid) for gid in genome_ids))
        nn_rows.append(dict(biosample=nid,
                             nearest_genome_1=dists[0][1], nearest_genome_1_distance=round(dists[0][0], 6),
                             nearest_genome_2=dists[1][1], nearest_genome_2_distance=round(dists[1][0], 6),
                             nearest_genome_3=dists[2][1], nearest_genome_3_distance=round(dists[2][0], 6),
                             median_distance_to_all_32_genomes=round(sorted(1.0 - jaccard(sketches[nid], sketches[g]) for g in genome_ids)[16], 6)))
    with open('results/g7_followup/g7f_structure_nearest_neighbors.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, list(nn_rows[0].keys()), delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(nn_rows)

    with open('results/g7_followup/g7f_structure_timing.json', 'w') as fh:
        json.dump({'kmer_size': KMER_SIZE, 'sketch_size': SKETCH_SIZE, 'hash': 'zlib.crc32',
                   'method': 'SCREENING_STRUCTURE_ONLY_read_based_bottom_minhash_bounded_subsample',
                   'read_pairs_cap_per_sample': READ_PAIRS_CAP, 'seconds_per_sample': timing}, fh, indent=2)
    print('done', len(all_ids), 'total points (32 genomes + 10 read-based samples)')


if __name__ == '__main__':
    main()
