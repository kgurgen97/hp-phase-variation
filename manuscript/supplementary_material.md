# Supplementary material

## Supplementary Methods

### S1. Reproducibility boundary

The publication package separates tracked manuscript-facing outputs from raw public sequencing data and large regenerable intermediates. Frozen generating-source exports are retained under `analysis/frozen_source/`, while normalized publication outputs are retained under `metadata/` and `results/`. The exact historical development environment was only partially recorded; the publication package therefore declares a newly tested Python 3.11 baseline rather than retrospectively reconstructing undocumented historical package versions.

The authoritative Phase-1 call matrix is not duplicated in the clean publication repository. Its development-archive path is `results/g5/g5_phase1_calls.tsv`, Git blob SHA is `58486570113e1cde6d0668db1b8d2a6fe38b1915`, and file size is 17,870,238 bytes. The clean repository records this identity together with the regeneration script and raw-data manifest in `results/calls/phase1_calls.MANIFEST.tsv`.

### S2. Reference-panel construction and repeat discovery

Repeat discovery used a frozen eight-genome *H. pylori* RefSeq panel. Reference FASTA and GFF accessions, checksums, and source URLs are tracked in `metadata/loci/reference_panel.tsv`. The publication analysis treats this panel as a pragmatic multi-reference discovery set rather than as exhaustive representation of global *H. pylori* diversity.

Frozen discovery rules considered 1–6 nt motifs. Pure homopolymers required at least 8 bases; dinucleotide repeats required at least 5 copies; 3–4 nt motifs required at least 4 copies; and 5–6 nt motifs required at least 3 copies with total tract length of at least 12 bp. Motifs were normalized to the shortest canonical period across rotations and strands. Loci longer than 100 bp were recorded but not caller-eligible. Unique anchoring required the immediate 25-bp left and right flanks to be unique in the corresponding genome and not themselves contain a qualifying repeat tract.

### S3. Repeat-locus evidence terminology

Repeat-locus membership, phase-variation evidence, and technical eligibility are separate variables. Member-level evidence includes REPEAT_ONLY and STRONG_PV_CANDIDATE; literature-anchored loci can separately carry KNOWN_PV evidence. The frozen technical classes PRIMARY_TECHNICAL, PROVISIONAL_TECHNICAL, and PRIMARY_INELIGIBLE are measurement/analysis classes and do not imply experimentally demonstrated ON/OFF switching.

### S4. Frozen caller gates and validation denominators

Read extraction used 25-bp flanking anchors and a 12-mer seed index. Up to two mismatches per anchor were allowed during extraction. At call stage, supporting fragments required mean anchor quality ≥20 and at most three total anchor mismatches. A locus required at least eight filtered spanning fragments, a dominant-allele fraction ≥0.40, a top-two-allele fraction ≥0.80, and an off-unit fraction ≤0.30 to remain callable. Confidence classification incorporated depth, strand balance, and the Wilson lower confidence bound on dominant fraction. The frozen noise model included motif-class-specific stutter estimates derived during development and a twofold safety factor. These are measurement rules, not disease-derived thresholds.



The synthetic benchmark used a deterministic disease-blind locus subset stratified by motif class and CDS/non-CDS context, with at most one simulated catalogue member per locus. Reads were simulated as SE150 or PE250 at 10×, 30×, 100×, and 300× nominal depth and at 100:0, 90:10, 70:30, and 50:50 allele-mixture ratios. Empirical stutter parameters estimated from invariant callable loci in PRJNA622860 were used for the authoritative EMPIRICAL mode. Extraction and calling used the full catalogue index. Synthetic performance is reported with explicit denominator nesting. The empirical benchmark contained 2,720 comparisons, of which 1,780 were callable. Exact dominant-repeat-length concordance is therefore reported as 1,757/1,780 among callable comparisons, while false-confident error is reported as 15/1,595 among HIGH/MEDIUM-confidence callable comparisons.

An earlier eight-sample ONT Stage A comparison was preserved as INCONCLUSIVE_VALIDATION for exact repeat-length inference because homopolymer ambiguity and a directional one-base discrepancy could not be attributed uniquely to ONT, the short-read caller, or shared locus boundaries; the G5 caller was not retuned. The high-accuracy rescue therefore used five isolates from five distinct patients selected by a frozen disease-blind HiFi-yield rule from CRA026546, with same-BioSample PacBio Sequel II CCS and MGISEQ-2000 paired short reads. Authors' short-read-polished assemblies were not used as truth.

For HiFi extraction, the tracked executed generator uses 25-bp catalogue anchors, exact 16-mer seed hits sampled across the read, Hamming verification with at most two mismatches per anchor, both orientations, and tract spans up to 150 bp. HIGH_CONFIDENCE_HIFI_TRUTH requires at least 20 spanning reads, a dominant fraction ≥0.90, unit alignment for motif length ≥2, and at least two dominant-allele reads per strand. Fewer than 10 spanning reads is UNCALLABLE_HIFI; other cases with at least 10 spanning reads that fail high-confidence or mixture criteria remain AMBIGUOUS_HIFI. The full isolate-locus truth universe comprised 5,980 pairs: 1,034 HIGH_CONFIDENCE_HIFI_TRUTH, 1,453 AMBIGUOUS_HIFI, and 3,493 UNCALLABLE_HIFI. Short-read exact concordance was assessed in the high-confidence comparison universe (911/913). The non-reference subset contained 84 comparisons, of which 82 were exact. Motif strata with zero observed non-reference truth are not claimed to be validated for non-reference allele recovery.

### S5. Pilot endpoint gate

Primary disease analysis was restricted to the 261 PRIMARY_TECHNICAL mono-A/T loci. Before labels were used for testing, the frozen per-contrast locus screen required: (i) overall callability fraction ≥0.80 across the contrast sample set; (ii) at least two distinct observed dominant repeat lengths among callable samples; and (iii) a minor-allele count ≥2. Failing loci were excluded without imputation or recoding. PROVISIONAL_TECHNICAL loci were processed separately and never merged into the primary multiple-testing family.

Pairwise distance was the number of discordant dominant repeat-length alleles divided by the number of jointly callable screened PRIMARY loci. The strict cohort-level multivariate endpoint required at least 10 jointly callable loci for every sample pair entering the distance analysis. Group separation statistics used exhaustive enumeration of all label allocations preserving observed group sizes. R² uncertainty used 2,000 stratified within-group bootstrap replicates with frozen Python random seed 88022026. PRJNA360417 had four primary loci passing the frozen cohort-level screen, with pairwise coverage 2–4 loci. PRJNA678459 had six primary loci passing the screen (plus one provisional locus), with pairwise coverage 4–6 loci. Both strict endpoints were therefore NOT_COMPUTABLE. Below-floor permutation statistics are retained only as diagnostics.

PRJNA1103397 was GC-only in the frozen pilot package and had zero primary loci passing the cohort-level screen; no disease contrast was performed.

PRJNA678459 used the already-frozen disease-blind G7-follow-up caller output, generated with extraction capped at 500,000 read fragments per sample. This bounded execution was retained unchanged in G8 and is not described as a full-depth rerun.

### S6. Structure diagnostic

PRJNA678459 received a disease-blind read-based k-mer/MinHash screening structure analysis. The frozen implementation used canonical k=21, a bottom-1,000 sketch with deterministic zlib.crc32 hashing, and the first 30,000 read pairs per sample. The 10 read-based sample sketches were positioned relative to the frozen 32-genome G7 screening panel. Disease labels were not read by the structure-generation script and were overlaid only in G8. The frozen result was pseudo-F 1.0506, R² 0.1161, exact p 0.031746 with 252 permutations. This is reported only as a confounding warning. It is not an ancestry-adjusted disease model and is not evidence for a repeat-length disease association. PRJNA360417 remains NOT_PERFORMED for structure analysis.


### S7. HiFi extractor implementation provenance note

The frozen high-accuracy truth rules file contains one stale implementation-level description in rule H05: it states “exact 11-mer seeding plus banded semi-global edit distance”. The tracked generator that produced the frozen HiFi evidence and was committed with the G6 rescue instead implements exact 16-mer seeding followed by direct 25-bp anchor Hamming verification with at most two mismatches. The G6 rescue report explicitly describes the executed extractor as an own implementation using exact k-mer seeding with Hamming anchors. We therefore report the tracked generator as the executed implementation and retain the H05 text unchanged as historical provenance. The truth-status thresholds H06–H13, sample selection, frozen short-read caller, comparison universe, and published frozen result tables are unaffected by this documentation discrepancy. No re-extraction or rule change was performed during manuscript preparation.

## Supplementary Tables

### Table S1. Frozen eight-genome reference panel

| Order | Short code | Strain | RefSeq assembly | Assembly |
|---:|---|---|---|---|
| 1 | 26695 | 26695 | GCF_000008525.1 | ASM852v1 |
| 2 | J99 | J99 | GCF_000008785.1 | ASM878v1 |
| 3 | HPAG1 | HPAG1 | GCF_000013245.1 | ASM1324v1 |
| 4 | G27 | G27 | GCF_000021165.1 | ASM2116v1 |
| 5 | P12 | P12 | GCF_000021465.1 | ASM2146v1 |
| 6 | B8 | B8 | GCF_000196755.1 | ASM19675v1 |
| 7 | Shi470 | Shi470 | GCF_000020245.1 | ASM2024v1 |
| 8 | India7 | India7 | GCF_000185185.1 | ASM18518v1 |

Full checksums and source URLs are retained in `metadata/loci/reference_panel.tsv`.

### Table S2. Frozen catalogue and technical-universe counts

| Metric | Count |
|---|---:|
| Repeat loci | 1,405 |
| Catalogue members | 4,009 |
| Caller-eligible loci | 1,196 |
| BOUNDARY_STABLE | 400 |
| BOUNDARY_AMBIGUOUS | 48 |
| ORTHOLOGY_AMBIGUOUS | 497 |
| INSUFFICIENT_TO_RESOLVE | 460 |
| PRIMARY_TECHNICAL | 261 |
| PROVISIONAL_TECHNICAL | 39 |
| PRIMARY_INELIGIBLE | 1,105 |
| REPEAT_ONLY within PRIMARY | 243 |
| STRONG_PV_CANDIDATE within PRIMARY | 15 |
| Mixed evidence within PRIMARY | 3 |
| KNOWN_PV within PRIMARY | 0 |

### Table S3. Technical validation summary

| Validation quantity | Numerator | Denominator | Fraction |
|---|---:|---:|---:|
| Synthetic callable | 1,780 | 2,720 | 65.4% |
| Synthetic exact among callable | 1,757 | 1,780 | 98.7% |
| Synthetic confident among callable | 1,595 | 1,780 | 89.6% |
| Synthetic false-confident wrong | 15 | 1,595 | 0.94% |
| HiFi exact, high-confidence comparisons | 911 | 913 | 99.8% |
| HiFi exact, non-reference comparisons | 82 | 84 | 97.6% |

### Table S4. Non-reference HiFi validation by motif stratum

| Stratum | Non-reference comparisons | Exact | Interpretation |
|---|---:|---:|---|
| MONO_AT | 80 | 78 | Observed non-reference truth |
| MONO_GC | 3 | 3 | Observed non-reference truth |
| DINUCLEOTIDE | 0 | 0 | Non-reference alleles not validated |
| TRINUCLEOTIDE | 1 | 1 | Sparse non-reference truth |
| MOTIF_GE4 | 0 | 0 | Non-reference alleles not validated |
| KNOWN_PV | 0 | 0 | Non-reference alleles not validated |

### Table S5. Pilot disease-cohort endpoint status

| Cohort | Frozen contrast | Primary loci passing screen | Provisional loci passing screen | Pairwise coverage | Strict endpoint | Diagnostic exact p |
|---|---|---:|---:|---|---|---:|
| PRJNA360417 | NAG vs IM (5 vs 6) | 4 | 0 | 2–4 | NOT_COMPUTABLE | 0.770563 |
| PRJNA678459 | AG vs GC (5 vs 5), exploratory | 6 | 1 | 4–6 | NOT_COMPUTABLE | 0.142857 |
| PRJNA1103397 | GC-only descriptive | 0 | NA | NA | No contrast | NA |

Diagnostic p-values are below the pre-specified callability floor and are not primary inferential endpoints.

### Table S6. Reproduction map

| Result family | Frozen source | Publication output | Boundary |
|---|---|---|---|
| Catalogue construction | `analysis/frozen_source/g5/` | `metadata/loci/repeat_catalogue.tsv` | Exact reference assemblies and historical inputs required for full rebuild |
| Synthetic benchmark | G5 benchmark scripts + caller modules | `results/validation/` | Aggregate authoritative counts tracked |
| Phase-1 callability | G5 calling/callability scripts | `results/callability/` | Public reads required; full call matrix manifest-backed |
| Boundary/orthology audit | G6 audit scripts | `metadata/validation/catalogue_boundary_audit.tsv` | Validation reads needed for independent read evidence |
| HiFi validation | G6 HiFi scripts | `results/validation/` | Public HiFi/MGI reads required; tracked extractor is the executed implementation; H05 prose discrepancy documented in S7 |
| Technical universe | G7 universe script | `metadata/loci/analysis_universe.tsv` | Upstream curated closure-table provenance disclosed |
| Confounding audit | G7 confounding script | `results/confounding/` | Full harmonization history retained in development archive |
| PRJNA678459 bounded caller input | G7-follow-up caller/consolidation scripts | `results/calls/prjna678459_locus_calls.tsv` | Public reads required; frozen extraction cap 500,000 fragments/sample |
| Pilot disease analysis | G8 pipeline | `results/pilot_disease/` | PRJNA360417 requires Phase-1 matrix regeneration; PRJNA678459 consumes the frozen bounded G7-follow-up calls |
| PRJNA678459 structure generation | G7 k-mer/MinHash + G7-follow-up structure generator | `results/structure/` | Public reads/assemblies required; frozen read cap 30,000 pairs/sample; screening only |
| Structure diagnostic overlay | G8 structure script | `results/pilot_disease/structure_overlay_summary.tsv` | PRJNA678459 descriptive overlay only; PRJNA360417 NOT_PERFORMED |
