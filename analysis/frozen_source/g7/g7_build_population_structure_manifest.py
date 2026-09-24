#!/usr/bin/env python3
"""G7 section 3: build the disease-blind population-structure manifest from
small public assemblies already present locally (no raw-read download in
this gate; no new large download performed by this script).

Two sources, both already downloaded to data/reference/ (gitignored, small
assembly FASTA only):
  1. data/reference/g7_panel/ - the same 8-genome reference panel used to
     build the G5 repeat catalogue (metadata/g5/g5_reference_panel.tsv);
     used only as diverse reference/context isolates for lineage anchoring
     (role = REFERENCE_LINEAGE_ANCHOR), never as disease-labelled samples.
  2. data/reference/g7_prjna521871/ - 24 assembly-only BioSamples from
     PRJNA521871 (Nicaragua; Correa AG/IM disease labels from the BioSample
     host_disease field; patient_uid confirmed via host_subject_id in G2/G3
     harmonization, metadata/g3/g3_candidate_biosamples.tsv). This is the
     only currently source-checked disease-labelled cohort with public
     assembly sequence (not raw reads), so it is the only disease-labelled
     role in this manifest; see reports/g7/G7_COHORT_RESCUE.md for why the
     five newly-prioritized cohorts could not be added here (none have
     public assemblies; all are read-only and raw-read download is not
     authorized in this gate).

No FASTQ/BAM/CRAM/SRA/fast5/pod5 download performed here.

Run: python3 scripts/g7/g7_build_population_structure_manifest.py
"""
import csv
import hashlib
import os

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PANEL_TSV = os.path.join(REPO, "metadata", "g5", "g5_reference_panel.tsv")
INVENTORY_TSV = os.path.join(REPO, "metadata", "raw", "dataset_inventory.tsv")
G3_BIOSAMPLES_TSV = os.path.join(REPO, "metadata", "g3", "g3_candidate_biosamples.tsv")
PANEL_DIR = os.path.join(REPO, "data", "reference", "g7_panel")
PRJNA521871_DIR = os.path.join(REPO, "data", "reference", "g7_prjna521871")
OUT_TSV = os.path.join(REPO, "metadata", "g7", "g7_population_structure_manifest.tsv")

FIELDS = [
    "sample_id", "biosample", "patient_uid", "bioproject", "assembly_accession",
    "disease_label", "disease_label_provenance", "country", "region",
    "platform", "biopsy_site", "study", "role", "local_fasta_path",
    "local_sha256", "source_url", "inclusion_exclusion_reason",
]


def sha256_of(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def build_panel_rows():
    with open(PANEL_TSV, newline="") as fh:
        panel = list(csv.DictReader(fh, delimiter="\t"))
    files = {f.split(".")[0]: f for f in os.listdir(PANEL_DIR)}
    rows = []
    for p in panel:
        acc = p["assembly_accession_version"].split(".")[0]
        fname = files.get(acc)
        local_path = os.path.join("data", "reference", "g7_panel", fname) if fname else "NA"
        sha = sha256_of(os.path.join(REPO, local_path)) if fname else "NA"
        rows.append({
            "sample_id": f"REF_{p['strain']}",
            "biosample": "NA",
            "patient_uid": "NOT_APPLICABLE_REFERENCE_STRAIN",
            "bioproject": "NA",
            "assembly_accession": p["assembly_accession_version"],
            "disease_label": "NOT_APPLICABLE_REFERENCE_STRAIN",
            "disease_label_provenance": "NOT_APPLICABLE",
            "country": "UNKNOWN_REFERENCE_STRAIN",
            "region": "UNKNOWN",
            "platform": "NA_ASSEMBLY_ONLY",
            "biopsy_site": "UNKNOWN",
            "study": "G5_REFERENCE_PANEL",
            "role": "REFERENCE_LINEAGE_ANCHOR",
            "local_fasta_path": local_path,
            "local_sha256": sha,
            "source_url": p["fasta_url"],
            "inclusion_exclusion_reason": "INCLUDED: diverse global reference strain used for G5 catalogue; lineage anchoring only, disease-blind",
        })
    return rows


def build_prjna521871_rows():
    with open(INVENTORY_TSV, newline="") as fh:
        inv = [r for r in csv.DictReader(fh, delimiter="\t") if r["bioproject"] == "PRJNA521871"]
    with open(G3_BIOSAMPLES_TSV, newline="") as fh:
        g3 = {r["biosample"]: r for r in csv.DictReader(fh, delimiter="\t") if r["bioproject"] == "PRJNA521871"}
    files = {}
    for f in os.listdir(PRJNA521871_DIR):
        acc = f.split(".")[0]
        files[acc] = f
    rows = []
    for r in inv:
        acc = r["assembly_accession"]
        fname = files.get(acc)
        local_path = os.path.join("data", "reference", "g7_prjna521871", fname) if fname else "NA"
        sha = sha256_of(os.path.join(REPO, local_path)) if fname else "NA"
        g3row = g3.get(r["biosample"], {})
        disease = g3row.get("disease_harmonized", "UNKNOWN")
        reason = (
            "INCLUDED: disease-labelled (BioSample host_disease field, EXPLICIT for AG/IM), "
            "patient_uid CONFIRMED via host_subject_id, public assembly available; "
            "publication not found (REPOSITORY_ONLY) so disease label rests on the BioSample "
            "field alone, not a peer-reviewed sample-level table -- used for population-structure "
            "screening only, not primary disease inference"
            if fname else
            "EXCLUDED: no local assembly file matched for this accession"
        )
        if disease == "UNKNOWN":
            reason = (
                "INCLUDED as UNKNOWN-disease context (host_disease = 'Gastritis' without qualifier, "
                "not mapped to a Correa stage per G2 rules); structure-only, no disease label used"
                if fname else "EXCLUDED: no local assembly file matched for this accession"
            )
        rows.append({
            "sample_id": r["biosample"],
            "biosample": r["biosample"],
            "patient_uid": g3row.get("patient_uid", "NA"),
            "bioproject": "PRJNA521871",
            "assembly_accession": acc,
            "disease_label": disease,
            "disease_label_provenance": g3row.get("disease_evidence_source", "BioSample.host_disease"),
            "country": "Nicaragua",
            "region": r.get("geographic_region", "UNKNOWN"),
            "platform": "NA_ASSEMBLY_ONLY",
            "biopsy_site": g3row.get("biopsy_site_harmonized", "UNKNOWN"),
            "study": "PRJNA521871",
            "role": "DISEASE_LABELED_CONTEXT_ASSEMBLY" if disease != "UNKNOWN" else "UNKNOWN_DISEASE_CONTEXT_ASSEMBLY",
            "local_fasta_path": local_path,
            "local_sha256": sha,
            "source_url": f"https://www.ebi.ac.uk/ena/browser/view/{acc}",
            "inclusion_exclusion_reason": reason,
        })
    return rows


def main():
    rows = build_panel_rows() + build_prjna521871_rows()
    os.makedirs(os.path.dirname(OUT_TSV), exist_ok=True)
    with open(OUT_TSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    n_with_file = sum(1 for r in rows if r["local_fasta_path"] != "NA")
    print("rows:", len(rows), "with local fasta:", n_with_file)


if __name__ == "__main__":
    main()
