#!/usr/bin/env python3
"""G7 follow-up: build a QC manifest for the 10 verified local PRJNA678459 runs, same schema as the frozen
metadata/g4/g4_phase1_download_manifest.tsv, disease-blind (no disease/stage column; see g7f_scope_rules.tsv F06).
Consumes only the already-verified metadata/g7_followup/g7f_frozen_membership_manifest.tsv and
metadata/harmonized/run_linkage.tsv (patient_mapping_status only, not disease)."""
import csv, os

RAW = 'data/raw/g7_followup/PRJNA678459'
link = {r['run_accession']: r for r in csv.DictReader(open('metadata/harmonized/run_linkage.tsv'), delimiter='\t') if r['bioproject'] == 'PRJNA678459'}
mem = {r['run_accession']: r for r in csv.DictReader(open('metadata/g7_followup/g7f_frozen_membership_manifest.tsv'), delimiter='\t')}

cols = ['bioproject', 'biosample', 'run_accession', 'cohort_role', 'patient_unit_or_basis', 'library_layout', 'instrument',
        'fastq_url_1', 'fastq_url_2', 'fastq_url_unpaired', 'expected_bytes', 'expected_bytes_1', 'expected_bytes_2', 'expected_bytes_unpaired',
        'expected_md5_1', 'expected_md5_2', 'expected_md5_unpaired', 'ena_read_count', 'ena_base_count',
        'local_path_1', 'local_path_2', 'local_path_unpaired', 'download_status', 'checksum_status', 'actual_bytes']
rows = []
for run, r in sorted(link.items()):
    m = mem[run]
    d = os.path.join(RAW, run)
    f1 = os.path.join(d, run + '_1.fastq.gz')
    f2 = os.path.join(d, run + '_2.fastq.gz')
    rows.append(dict(bioproject='PRJNA678459', biosample=r['biosample'], run_accession=run,
                      cohort_role='G7_FOLLOWUP_TARGETED_RAW_READ_PROSPECTIVE_REPLICATION_CANDIDATE',
                      patient_unit_or_basis='CONFIRMED_EXACT_MAPPING_REPOSITORY_LEVEL_ONLY_NO_PUBLICATION (see metadata/g7/g7_disease_cohort_verification.tsv)',
                      library_layout='PAIRED', instrument=r['instrument'],
                      fastq_url_1='NA_LOCAL_ONLY_NO_NETWORK_THIS_RUN', fastq_url_2='NA_LOCAL_ONLY_NO_NETWORK_THIS_RUN', fastq_url_unpaired='NA',
                      expected_bytes=m['expected_bytes_est'], expected_bytes_1='NA', expected_bytes_2='NA', expected_bytes_unpaired='NA',
                      expected_md5_1='NOT_AVAILABLE_NO_NETWORK_ACCESS_THIS_RUN', expected_md5_2='NOT_AVAILABLE_NO_NETWORK_ACCESS_THIS_RUN', expected_md5_unpaired='NA',
                      ena_read_count='NOT_AVAILABLE_NO_NETWORK_ACCESS_THIS_RUN', ena_base_count='NOT_AVAILABLE_NO_NETWORK_ACCESS_THIS_RUN',
                      local_path_1=f1, local_path_2=f2, local_path_unpaired='NA',
                      download_status='PRE_EXISTING_LOCAL_REUSED_NOT_DOWNLOADED_THIS_RUN',
                      checksum_status=('PASS_SIZE_AND_GZIP_ONLY' if m['run_status'] == 'PASS_SIZE_AND_GZIP_ONLY' else 'FAIL'),
                      actual_bytes=m['actual_bytes']))

with open('metadata/g7_followup/g7f_qc_manifest.tsv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, cols, delimiter='\t', lineterminator='\n')
    w.writeheader()
    w.writerows(rows)
print('wrote', len(rows), 'rows')
