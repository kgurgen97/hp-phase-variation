# A reproducible short-read framework for profiling repeat-length variation at candidate phase-variable loci in *Helicobacter pylori*

## Abstract

Repeat-length variation can provide a reversible source of within-species phenotypic diversity in *Helicobacter pylori*, but short-read analysis is technically difficult: repetitive sequence, uncertain locus boundaries, orthology ambiguity, and sequencing noise complicate allele measurement, while population structure and cohort design can confound downstream disease association. We developed a conservative, reproducible short-read framework for cataloguing and calling repeat-length alleles at candidate phase-variable loci and evaluated its technical performance before applying it to public gastric-disease cohorts. An eight-genome reference panel yielded 1,405 repeat loci represented by 4,009 catalogue members, of which 1,196 loci had at least one caller-eligible member. A boundary and orthology audit classified 400 loci as boundary-stable, 48 as boundary-ambiguous, 497 as orthology-ambiguous, and 460 as insufficient to resolve. The frozen technical analysis universe comprised 261 PRIMARY_TECHNICAL loci, 39 PROVISIONAL_TECHNICAL loci, and 1,105 PRIMARY_INELIGIBLE loci. Importantly, technical eligibility was kept distinct from biological phase-variation evidence: the PRIMARY_TECHNICAL set contained 243 repeat-only loci, 15 strong phase-variation candidates, three loci with mixed member-level evidence, and no literature-anchored KNOWN_PV loci. In an empirical synthetic benchmark, 1,780 of 2,720 comparisons were callable; 1,757/1,780 callable comparisons matched the expected dominant repeat-length allele (98.7%), while 15/1,595 high- or medium-confidence calls were confidently wrong (0.94%). Against high-confidence HiFi truth, 911/913 short-read comparisons were exact (99.8%); among non-reference alleles, 82/84 were exact. Public-cohort evaluation showed strong coupling of study and geography together with study/platform-associated callability differences. Under pre-specified pilot rules, the two disease contrasts with sufficient clinical labels failed the minimum joint-callability floor and therefore had non-computable strict endpoints. The framework therefore supports conservative repeat-length profiling from short reads, while the available public data do not support a definitive test of repeat-state association across gastric carcinogenesis.

## Impact statement

Repeat tracts are biologically important in *Helicobacter pylori*, but their instability also makes them difficult to measure and interpret from short-read sequencing. This study separates repeat discovery, evidence for phase variation, technical callability, and downstream disease-analysis eligibility rather than treating them as interchangeable. The resulting framework shows high exact concordance in callable synthetic and HiFi-validated comparisons while explicitly retaining uncallable loci, orthology uncertainty, and cohort confounding. Applying the frozen framework to public gastric-disease datasets also provides an informative identifiability result: the available cohorts do not support computation of the originally intended strict disease-association endpoint under the pre-specified rules. The work therefore provides a reproducible technical foundation and defines the data requirements for a future disease-focused discovery-and-replication study.

## Data summary

No new patient recruitment or primary sequencing was performed for this study. All raw sequencing data were previously public. The principal data sources were NCBI BioProjects PRJNA360417 (Colombian gastric-lesion isolates) [15], PRJNA622860 (within-host *H. pylori* diversification; used for caller/noise development) [16], PRJNA816422 (paired short-read/ONT validation data) [17], PRJNA1103397 (Chinese gastric-cancer/gastritis isolates) [18], and PRJNA678459 (10 Gansu, China isolates used only in the bounded exploratory follow-up), together with CNCB-NGDC BioProject PRJCA041148 / Genome Sequence Archive CRA026546 for the PacBio HiFi/MGISEQ rescue set [19]. For PRJNA678459, the NCBI BioProject record was the authoritative repository source; no linked source publication was identified in our publication audit.

Individual BioSample and run accessions are enumerated in the tracked manifests under `metadata/cohorts/`, `metadata/validation/`, and the frozen G7-follow-up manifests. The publication repository contains the caller, frozen rules, catalogue and eligibility tables, derived manuscript-facing outputs, numerical audits, and figure-rendering code. The large Phase-1 call matrix is not duplicated in the clean repository; its authoritative development-archive path, Git blob SHA, byte size, raw-data manifest, and regeneration route are recorded in `results/calls/phase1_calls.MANIFEST.tsv`. A persistent public archive/DOI for the publication repository should be created before final submission or acceptance, according to journal data-availability requirements.

## Keywords

*Helicobacter pylori*; phase variation; tandem repeats; short-read sequencing; long-read validation; population structure; gastric carcinogenesis; reproducible genomics

## Data Summary

All manuscript-facing code, configuration, metadata, validation summaries, pilot outputs, numerical audits, and figure-rendering scripts are contained in the publication repository. Public sequencing accessions used by the frozen workflow are listed in the tracked manifests under `metadata/cohorts/` and `metadata/validation/`. Raw sequencing data are not redistributed.

The large Phase-1 call matrix is not duplicated in the clean publication repository. Its authoritative development-archive location, Git blob SHA (`58486570113e1cde6d0668db1b8d2a6fe38b1915`), size (17,870,238 bytes), reconstruction script, and public raw-data manifest are recorded in `results/calls/phase1_calls.MANIFEST.tsv`. Manuscript-facing derived values are stored in `results/manuscript/manuscript_numbers.tsv` and checked automatically against tracked frozen outputs.

## Impact Statement

Repeat-mediated phase variation is biologically important in *Helicobacter pylori*, but repetitive sequence and extreme bacterial population structure make short-read measurements unusually vulnerable to technical and epidemiological confounding. This study provides a reproducible framework that separates repeat-locus discovery, biological evidence for phase variation, technical callability, long-read validation, and disease-analysis eligibility. The caller showed high exact concordance in synthetic and HiFi-supported comparisons, including non-reference alleles, while the broader analysis exposed substantial boundary, orthology, and callability limitations. Applying frozen rules to public gastric-disease cohorts showed that the available data could not support the intended strict disease-association endpoint. The work therefore supplies a validated technical resource and, equally importantly, a transparent definition of when public short-read data are insufficient for biological inference.

## Introduction

*Helicobacter pylori* is a genetically diverse gastric bacterium in which repeat-length changes can alter coding or regulatory sequence and, for experimentally supported loci, generate reversible functional states [1–5]. This biology makes short tandem and homopolymeric repeat loci attractive candidates for studying rapid adaptation within the gastric niche. It also creates a measurement problem: the same sequence properties that make repeat tracts mutable make them difficult to infer reliably from short reads, a challenge also recognized in broader STR genotyping work [6,7].

Several distinct sources of error can mimic or obscure true repeat-length variation. Short-read alignment and local sequence context can distort repeat counts; orthologous loci can be difficult to define across diverse *H. pylori* backgrounds; reference-centred catalogues can fragment or merge biologically distinct loci; and sequencing platform or study design can become entangled with apparent callability. These issues are particularly important for disease-association analyses because *H. pylori* has pronounced phylogeographic population structure, documented from early multilocus studies through recent worldwide whole-genome collections [8–11]. A technically plausible allele difference is therefore not automatically a reliable phase-variable state, and a cohort-level difference is not automatically a disease-associated effect.

We designed the present study around that distinction. The primary objective was not to claim that every repeat locus is an experimentally established phase-variable locus, nor to force a gastric-cancer association analysis from heterogeneous public datasets. Instead, we sought to build and validate a conservative short-read framework that separates: (i) repeat-locus discovery, (ii) biological evidence for phase variation, (iii) technical callability, and (iv) downstream disease-analysis eligibility.

The study was developed in staged gates against the biological background of the Correa gastric precancerous cascade [10–12]. We first constructed a multi-reference repeat catalogue, audited locus boundaries and orthology, and defined technical eligibility. We then evaluated the caller using empirical synthetic data and independent high-accuracy long-read truth. Finally, with all analytical rules frozen, we tested whether available public gastric-disease cohorts could support the originally intended disease-association analysis across disease states related to the Correa model of gastric carcinogenesis [12]. This last stage was treated as a stress test of identifiability rather than as a requirement to produce a positive association result.

## Methods

### Study design and reproducibility principles

The analysis was organized as a staged workflow in which technical definitions and thresholds were frozen before the final disease-association inspection. Original metadata values were preserved, harmonized values were kept traceable to their sources, ambiguous disease labels were not reassigned without explicit support, and multiple isolates from one patient were not treated as independent patients where patient identity was available or could be resolved.

The publication repository contains the frozen caller implementation, catalogue and eligibility tables, validation summaries, public-data manifests, pilot-analysis rules, manuscript-number audit, and publication-figure rendering code. Generating scripts from the development archive were exported verbatim for the catalogue, synthetic benchmark, Phase-1 callability, boundary/orthology audit, HiFi validation, technical locus universe, confounding audit, the bounded disease-blind G7 PRJNA678459 follow-up, pilot disease analysis, and structure diagnostic. Debug scripts and superseded internal orchestration were excluded from the publication export.

The historical development environment was only partially recorded. We therefore do not retrospectively claim an exact lockfile for every development gate. The publication package instead declares a newly tested Python 3.11 baseline and runs consistency checks on Python 3.10-3.12. From-raw regeneration additionally requires the public sequencing data and stage-specific external tools.

### Multi-reference repeat catalogue

Repeat discovery was performed across a frozen eight-genome *H. pylori* reference panel. Motifs of 1–6 nt were considered. Pure homopolymers required at least 8 bases; dinucleotide repeats required at least 5 copies; 3–4 nt motifs required at least 4 copies; and 5–6 nt motifs required at least 3 copies and a total tract length of at least 12 bp. Motifs were canonicalized across rotations and strands, and loci longer than 100 bp were retained in the catalogue but were not caller-eligible. Candidate repeat tracts were reconciled across references into locus-level groups using frozen flank- and annotation-based equivalence rules. The resulting catalogue records both a locus identifier and member-level reference representation.

The final catalogue contained 1,405 loci represented by 4,009 catalogue members. A locus was counted as caller-eligible when at least one of its catalogue members satisfied the caller-eligibility criteria. This produced 1,196 caller-eligible loci.

Biological evidence and technical eligibility were represented as separate fields. Catalogue members could be annotated as REPEAT_ONLY or STRONG_PV_CANDIDATE, and literature-anchored loci could separately carry KNOWN_PV evidence. These annotations were not used as synonyms for the final technical analysis class.

### Boundary and orthology audit

Because repeat tracts can be difficult to reconcile across divergent genomes, each locus underwent a frozen boundary/orthology audit. The audit assigned one of four states: BOUNDARY_STABLE, BOUNDARY_AMBIGUOUS, ORTHOLOGY_AMBIGUOUS, or INSUFFICIENT_TO_RESOLVE.

The audit was explicitly conservative. For multi-member loci, members were required to share a canonical motif and coarse CDS/non-CDS context, and tract-length differences for motifs of at least 2 nt had to be compatible with motif periodicity. Adjacent 25-bp flanks were compared in both relative orientations with small shifts; a locus required a connected member graph in which both flanks of each linking pair reached the frozen 0.72 identity threshold. Orthology ambiguity flags, product/neighbour-only grouping without flank support, or multiple members from the same reference genome triggered ORTHOLOGY_AMBIGUOUS; single-member loci lacked cross-member evidence and were classified INSUFFICIENT_TO_RESOLVE. A locus could therefore be technically measurable in one reference context yet remain unsuitable for cross-genome interpretation if repeat boundaries or orthology could not be resolved reliably.

### Short-read repeat-length caller

The caller estimates repeat-length alleles from fragments spanning the repeat tract between flanking sequence anchors. The frozen implementation uses 25-bp anchors, permits at most two mismatches per anchor during extraction and at most three total anchor mismatches at the call stage, requires mean anchor base quality of at least 20, and requires at least eight filtered spanning fragments for a callable locus. The two most frequent unit-aligned alleles must account for at least 80% of supporting reads and the dominant allele for at least 40%; excess off-unit signal or diffuse support renders the locus uncallable. Confidence classes additionally use read depth, strand balance, the Wilson lower bound on dominant-allele fraction, and a frozen stutter/noise model estimated during development. Real-sample outputs are interpreted as repeat-length alleles. They are not converted automatically into sample-specific ON/OFF functional states because functional consequences depend on locus-specific sequence context and biological validation.

### Synthetic benchmark

An empirical synthetic benchmark was used to evaluate exact dominant repeat-length recovery and confidence calibration under repeat-associated noise. Loci were selected deterministically and disease-blind across motif/context strata, with at most one catalogue member per locus; read generation was evaluated in SE150 and PE250 layouts, at nominal depths of 10×, 30×, 100×, and 300×, and at simulated allele-mixture ratios of 100:0, 90:10, 70:30, and 50:50. Empirical stutter parameters estimated from invariant callable loci in PRJNA622860 [16] were used for the authoritative empirical mode, while extraction and calling used the full catalogue index so cross-member ambiguity remained possible. The authoritative empirical benchmark comprised 2,720 comparisons. We report callability relative to all empirical comparisons, exact dominant-allele concordance among callable comparisons, confidence among callable comparisons, and the false-confident error fraction among high- or medium-confidence calls.

Thresholds were not altered after inspection of disease-association results.

### Independent HiFi validation

An initial Oxford Nanopore-based Stage A comparison used a disease-blind subset of PRJNA816422 [17] and was preserved as inconclusive for exact repeat-length validation because homopolymer ambiguity and a directional one-base discrepancy could not be attributed uniquely to the long-read platform, the short-read caller, or shared boundary definitions; the short-read caller was not retuned. We therefore used an independently selected high-accuracy rescue set from PRJCA041148 / CRA026546 [19], comprising five isolates from five distinct patients with same-BioSample PacBio Sequel II CCS/HiFi reads and MGISEQ-2000 paired short reads. PacBio circular-consensus/HiFi sequencing provides long reads with high per-read accuracy and is well suited to resolving repeat-containing sequence contexts [13].

The executed HiFi extractor located the same 25-bp catalogue anchors in raw CCS reads using exact 16-mer seeding followed by anchor Hamming verification with at most two mismatches, on both strands, and retained tract spans up to 150 bp with cross-locus arbitration. HIGH_CONFIDENCE_HIFI_TRUTH required at least 20 spanning HiFi reads, dominant-allele fraction at least 0.90, unit alignment for motifs of at least 2 nt, and at least two dominant-allele reads from each strand. Fewer than 10 spanning reads was UNCALLABLE_HIFI; remaining cases with at least 10 reads that did not meet high-confidence or pre-specified mixture criteria were AMBIGUOUS_HIFI. Short-read concordance was evaluated only in the callable/high-confidence comparison universe.

Because reference-length agreement alone can inflate apparent performance, a separate non-reference analysis quantified exact concordance when the HiFi-supported allele differed from the catalogue reference length. Motif-stratified summaries were retained to distinguish strata with observed non-reference truth from strata in which non-reference alleles were not represented.

### Frozen technical analysis universe

Technical locus classes were assigned after integrating catalogue properties, callability, boundary/orthology status, motif-class validation evidence, and frozen exclusion rules. PRIMARY_TECHNICAL required boundary stability, absence of the frozen study/platform-confounded callability label, caller eligibility/callability under the G5 rules, and an EMPIRICALLY_VALIDATED_EXACT motif-class interpretation; this yielded the 261 mono-A/T loci used in the primary disease analysis. Technically eligible loci from motif classes supported but limited by allele diversity were retained as PROVISIONAL_TECHNICAL, whereas all other loci were PRIMARY_INELIGIBLE. These classes were frozen before disease analysis.

PRIMARY_TECHNICAL should be interpreted as a technically selected repeat-length measurement universe. It is not an experimentally established set of phase-variable genes. PROVISIONAL_TECHNICAL loci were retained for sensitivity or descriptive use where appropriate but were not promoted into the primary universe.

### Public-cohort confounding audit

Before disease testing, public datasets with Correa-cascade-related metadata were audited for study, geography, sequencing platform, biopsy/disease context, and sample independence. The audit covered 45 studies and 793 Correa-labelled BioSamples in the frozen metadata table. Forty-one of 45 studies were single-country and four were multicountry.

Phase-1 callability was also compared across study/platform contexts. Loci labelled PLATFORM_CONFOUNDED_CALLABILITY are interpreted descriptively as study/platform-associated callability patterns; this label is not a causal estimate of platform effect.

### Pilot disease analysis

The disease application was deliberately restricted to cohorts that could support within-study contrasts under frozen metadata and independence rules. Primary disease analysis used only the 261 frozen PRIMARY_TECHNICAL mono-A/T loci; the 39 PROVISIONAL_TECHNICAL loci were kept in a separately labelled sensitivity layer and were never merged into the primary testing family. Three cohorts entered the final pilot package.

PRJNA360417 [15] provided a NAG-versus-IM contrast with 5 and 6 patients, respectively. PRJNA678459 provided an AG-versus-GC contrast with 5 and 5 repository-resolved independent sample units and was retained as a standalone exploratory, lower-provenance contrast. Its frozen caller input had been generated disease-blind from a bounded subsample of at most 500,000 read fragments per sample; the caller was not rerun or extended for the disease analysis. PRJNA1103397 [18] contributed six selected GC isolates; the source publication reported one selected isolate per patient, but explicit patient identifiers were unavailable in the frozen metadata. This cohort did not provide a usable within-cohort disease contrast in the frozen analysis.

Before disease labels were overlaid, each within-cohort contrast used the same frozen locus screen: overall callability of at least 0.80 across the contrast samples, at least two observed dominant repeat lengths among callable samples, and a minor-allele count of at least two. Failing loci were excluded without imputation or recoding. Pairwise repeat-allele distance was then defined as the fraction of jointly callable screened primary loci at which the two samples had different dominant repeat lengths. The pre-specified strict endpoint required at least 10 jointly callable loci for every sample pair retained in the distance matrix; otherwise the pair was undefined. Group separation used PERMANOVA statistics with exhaustive enumeration of all label allocations compatible with the observed group sizes, and R² uncertainty was summarized by 2,000 stratified bootstrap replicates with the frozen random seed. This floor and all screening rules were fixed before final disease-association inspection. When the floor was not met, the strict endpoint was declared NOT_COMPUTABLE. Below-floor permutation results, where available, were retained only as diagnostics and not promoted to inferential disease findings.

### Population-structure diagnostic

A disease-blind k-mer/MinHash distance analysis was used as a screening-level structure diagnostic for PRJNA678459; MinHash sketching provides an efficient approximation to genome-scale sequence distance [14]. The frozen implementation used canonical 21-mers, a bottom-1,000 sketch with deterministic zlib.crc32 hashing, and the first 30,000 read pairs per PRJNA678459 sample, positioning those samples relative to the existing 32-genome G7 screening panel. Disease labels were overlaid only after the distance matrix had been generated. This diagnostic was designed to identify whether disease labels aligned with detectable genomic structure; it was not treated as a phylogeny, an ancestry-adjusted association model, or a definitive confounding correction. PRJNA360417 did not receive a completed structure overlay and remains explicitly marked NOT_PERFORMED.

## Results

### A multi-reference catalogue produced a large repeat-locus universe but a much smaller technically interpretable subset

The frozen eight-genome panel produced 1,405 repeat loci represented by 4,009 catalogue members (Figure 1A). Of these, 1,196 loci had at least one caller-eligible member.

The boundary/orthology audit showed that cross-reference interpretation was a major constraint. Only 400 loci were classified BOUNDARY_STABLE, whereas 48 were BOUNDARY_AMBIGUOUS, 497 were ORTHOLOGY_AMBIGUOUS, and 460 were INSUFFICIENT_TO_RESOLVE (Figure 1B). Thus, apparent presence in a repeat catalogue did not imply that a locus could be transferred uncritically across strain backgrounds.

After all frozen technical filters, 261 loci entered PRIMARY_TECHNICAL, 39 entered PROVISIONAL_TECHNICAL, and 1,105 were PRIMARY_INELIGIBLE (Figure 1C).

The biological evidence composition of the PRIMARY_TECHNICAL set further illustrates why technical eligibility and phase-variation evidence must remain separate. Of the 261 primary loci, 243 were repeat-only, 15 were strong phase-variation candidates, and three contained mixed member-level evidence. No literature-anchored KNOWN_PV locus was present in PRIMARY_TECHNICAL (Figure 1D). The primary universe should therefore be read as a technically measurable repeat-locus panel, not as a validated phase-variable gene set.

### Synthetic benchmarking showed high exactness among callable comparisons but substantial uncallability

In the empirical synthetic benchmark, 1,780 of 2,720 comparisons were callable, leaving an uncallable fraction of 0.3456. Among callable comparisons, 1,757/1,780 recovered the exact dominant repeat-length allele, corresponding to 98.7% exact concordance (Figure 2A).

A total of 1,595 callable comparisons were assigned high or medium confidence. Fifteen of these 1,595 confident calls were wrong, corresponding to a false-confident error fraction of 0.94%. These results support a conservative interpretation: when the caller returns a technically supported allele, exact dominant-length recovery is high, but failure to call is common enough that missingness must be treated as a central property of the assay rather than ignored.

### Independent HiFi validation confirmed high concordance, including most observed non-reference alleles

Across the high-accuracy truth universe, 1,034 isolate-locus pairs were classified HIGH_CONFIDENCE_HIFI_TRUTH, 1,453 were AMBIGUOUS_HIFI, and 3,493 were UNCALLABLE_HIFI, for a total of 5,980 pairs (Figure 2B). The size of the ambiguous and uncallable fractions again shows that technical truth availability is itself selective.

Within the frozen high-confidence comparison universe, 911/913 short-read calls were exactly concordant with HiFi-supported repeat length (99.8%; Figure 2C). Importantly, 84 comparisons involved a HiFi-supported non-reference allele, and 82/84 of these were exact in the short-read call set (97.6%).

Non-reference validation coverage was uneven across motif strata (Figure 2D). Most non-reference observations occurred in mono-AT loci; some motif classes had no observed non-reference truth and therefore cannot be considered validated for non-reference allele recovery solely from the aggregate concordance value.

### Failure modes were dominated by locus-definition uncertainty and study/platform-associated callability

The catalogue audit demonstrated that boundary and orthology uncertainty affected the majority of repeat loci outside the stable subset (Figure 3A). This technical constraint is distinct from read-level calling performance: a read-level allele may be measured precisely while still being difficult to interpret across genomes if the underlying locus definition is unstable.

Phase-1 profiling identified a callable union of 909 loci across the frozen sample set. A total of 182 loci were labelled PLATFORM_CONFOUNDED_CALLABILITY in the frozen transfer analysis (Figure 3B). Because study, library construction, sequencing platform, geography, and biological cohort can be correlated in public datasets, these labels are treated as descriptive technical associations rather than causal platform effects.

The distribution of literature-anchored KNOWN_PV loci across technical classes also exposed an evidence gap (Figure 3C): none entered PRIMARY_TECHNICAL. This prevented us from presenting the primary technical universe as a validated set of phase-switching genes and motivated the terminology used throughout this manuscript.

Reproducibility boundaries are documented explicitly (Figure 3D). The exact historical development environment was incompletely recorded; the large Phase-1 call matrix is not duplicated in the clean publication repository but is identified by its authoritative archive path, Git blob SHA, byte size, raw-data manifest, and regeneration script; and selected G6 closure tables retain disclosed provenance limitations.

### Public gastric-disease cohorts did not support the intended strict disease endpoint

The confounding audit identified 45 Correa-labelled studies comprising 793 BioSamples. Forty-one studies were single-country and four were multicountry (Figure 4A), illustrating the strong coupling between study identity and geography in the public-data landscape.

The frozen pilot contrasts were small. PRJNA360417 contained 5 NAG and 6 IM patients; PRJNA678459 contained 5 AG and 5 GC patients; PRJNA1103397 contributed six selected GC isolates but no internal disease contrast (Figure 4B); explicit patient identifiers were unavailable in the frozen metadata.

Callability collapsed substantially when the strict locus screen was applied. In PRJNA360417, only 4/261 PRIMARY_TECHNICAL loci passed the cohort-level screen. In PRJNA678459, 6/261 primary loci and 1/39 provisional loci passed. PRJNA1103397 had 0/261 primary loci passing the screen (Figure 4C).

Neither disease contrast met the pre-specified requirement of at least 10 jointly callable loci per sample pair. Consequently, the strict endpoint was NOT_COMPUTABLE for both PRJNA360417 and PRJNA678459 (Figure 4D). The retained below-floor diagnostic exact permutation values were 0.770563 for PRJNA360417 and 0.142857 for PRJNA678459; these are diagnostics only and are not interpreted as formal disease-association tests.

For PRJNA678459, the screening structure analysis produced pseudo-F = 1.0506, R² = 0.1161, and exact p = 0.031746 from 252 label permutations. This result is treated as a confounding caution flag indicating alignment between disease labels and detectable genomic structure, not as evidence for or against a repeat-length disease association. The corresponding PRJNA360417 structure analysis remains NOT_PERFORMED.

## Discussion

This study establishes a conservative framework for measuring repeat-length variation in *H. pylori* short-read data while making explicit the technical and inferential boundaries that can otherwise be hidden by a single catalogue or association table.

The first major finding is methodological: catalogue construction is not equivalent to locus interpretability. Although the multi-reference panel yielded 1,405 repeat loci, only 400 were boundary-stable in the frozen audit, and the final primary technical universe contained 261 loci. Orthology ambiguity and unresolved boundaries therefore represent first-order constraints, not minor annotation details. For a highly diverse bacterial species, the locus-definition problem must be addressed alongside read-level allele calling.

The second finding is that the short-read caller performs strongly when high-confidence truth is available. Exact synthetic concordance was 1,757/1,780 among callable empirical comparisons, and HiFi concordance was 911/913 among high-confidence comparisons. The non-reference subset, 82/84 exact, is particularly important because reference-length agreement alone would provide a weak validation of true allele discrimination. At the same time, both synthetic and HiFi analyses showed substantial uncallability or truth ambiguity. A robust pipeline therefore needs to represent missingness explicitly rather than treating every catalogue locus as measurable in every sample.

The third finding concerns terminology and biological interpretation. The PRIMARY_TECHNICAL panel is not a synonym for experimentally established phase-variable loci. Most primary loci were repeat-only, only 15 were strong phase-variation candidates, three had mixed evidence, and none of the literature-anchored KNOWN_PV loci entered the primary technical set. For this reason, the sample-level observable throughout the analysis is a repeat-length allele. Functional ON/OFF state should only be assigned when locus-specific biology supports that mapping.

The public-cohort stress test also provides an important feasibility and identifiability result. The original biological question was whether reproducible phase-variable functional states associate with progression across gastric carcinogenesis after accounting for population structure and cohort effects. The available public data did not permit that question to be tested rigorously under the frozen rules. The limiting factors were not a lack of nominal p-values but insufficient jointly callable loci, small within-study disease contrasts, broader coupling between study, geography, and sequencing context, and—in PRJNA678459, the cohort with a completed structure diagnostic—alignment between disease labels and detectable genomic structure.

This distinction matters. A non-computable strict endpoint does not support a null biological conclusion, just as a below-threshold diagnostic p-value does not support a positive one. The central disease-association hypothesis therefore remains unresolved. The present contribution is instead a validated technical framework and an explicit demonstration of the data properties required for a future adequately powered discovery-and-replication study.

A future disease-focused study should be designed prospectively around patient-level independence, balanced disease stages within geography and lineage, consistent sequencing protocols, biopsy-site metadata, and sufficient depth to achieve a broad jointly callable locus set. Discovery and replication cohorts should be separated before final association testing, and ancestry adjustment should be incorporated directly into the inferential model rather than approximated by post hoc visualization.

### Limitations

Several limitations should be considered when using this framework.

First, technical validation is selective. High concordance is reported within callable or high-confidence truth subsets; it should not be extrapolated to loci or motif classes that remain uncallable or lack non-reference validation examples. The HiFi and short-read analyses also use the same catalogue anchor/boundary definition, so concordance does not independently validate that shared boundary definition; the separate disease-blind boundary/orthology audit addresses this limitation at the catalogue level.

Second, the frozen reference panel contains eight genomes and cannot represent the full global diversity of *H. pylori*. Orthology and boundary uncertainty observed here may therefore underestimate difficulties encountered in broader strain collections.

Third, the clean publication repository does not duplicate the 17.87 MB Phase-1 call matrix. The authoritative archive identity is preserved by repository, path, Git blob SHA (`58486570113e1cde6d0668db1b8d2a6fe38b1915`), byte size (17,870,238), raw-data manifest, and regeneration route. Manuscript-facing derived summaries are tracked and machine-audited.

Fourth, the exact historical development environment was not fully locked. The publication package provides a newly tested baseline rather than reconstructing undocumented historical package versions.

Fifth, the pilot disease cohorts are too small and technically sparse for definitive association testing. In addition, the frozen PRJNA678459 caller and structure inputs were intentionally bounded to 500,000 read fragments and 30,000 read pairs per sample, respectively, during the disease-blind G7 follow-up; they should not be represented as full-depth analyses. The diagnostic permutation statistics reported here must not be treated as primary inferential endpoints.

## Conclusions

A conservative short-read workflow can recover repeat-length alleles in *H. pylori* with high exact concordance when loci are technically callable and high-confidence truth is available. The more difficult problem is defining which loci are comparable across diverse genomes and which public cohorts can support unbiased downstream inference. By separating repeat discovery, phase-variation evidence, technical eligibility, validation, confounding assessment, and disease-analysis identifiability, this framework provides a reproducible foundation for future studies without overstating what current public data can establish.

The central hypothesis that phase-variable functional states reproducibly associate with progression across the Correa gastric carcinogenesis cascade remains unresolved. Testing it will require purpose-built discovery and independent replication datasets with adequate patient-level metadata, genomic structure control, and repeat-locus callability.

## Data and code availability

The publication repository contains the caller implementation, frozen configuration files, multi-reference repeat catalogue, technical locus universe, validation summaries, public-data manifests, confounding outputs, pilot disease outputs, generating-source export, manuscript-number audit, and publication-figure rendering code.

Raw sequencing data are not redistributed. Public accession manifests are retained for from-raw reconstruction. The large Phase-1 call matrix is retained in the authoritative development archive and is identified in the clean publication repository by exact archive path, Git blob SHA, file size, regeneration script, and raw-data manifest.

## Reproducibility statement

Core result counts and statistics used in the Abstract, Results, and main figures are recorded in `results/manuscript/manuscript_numbers.tsv`. The script `analysis/check_manuscript_numbers.py` verifies these values against tracked frozen artifacts; supplementary stratified values trace directly to the cited frozen TSV/JSON outputs. Publication figures are rendered from the frozen numerical table and tracked categorical outputs and do not recompute biological endpoints. The executed generating scripts are treated as the source of truth for implementation details; a documented discrepancy between the frozen H05 prose description and the tracked HiFi extractor implementation is retained explicitly in the Supplementary Methods rather than silently rewriting the frozen rules.

## Figure legends

### Figure 1. Repeat-locus catalogue and frozen technical analysis universe

(A) Catalogue scale: repeat loci, catalogue members, and caller-eligible loci. (B) Boundary/orthology audit states. (C) Final technical classes. (D) Evidence composition within PRIMARY_TECHNICAL. Technical eligibility is shown separately from phase-variation evidence.

### Figure 2. Technical validation of repeat-length calling

(A) Empirical synthetic benchmark, with numerators and denominators for callability, exact dominant-allele recovery, confidence, and false-confident error. (B) HiFi truth availability. (C) Exact short-read versus HiFi concordance overall and among non-reference alleles. (D) Non-reference validation coverage by motif stratum.

### Figure 3. Technical failure modes and interpretation boundary

(A) Catalogue-resolution states. (B) Study/platform-associated callability classes. (C) Distribution of literature-anchored KNOWN_PV loci across technical classes. (D) Reproducibility boundaries, including the historical environment, large Phase-1 matrix, and curated closure-table provenance.

### Figure 4. Public-cohort stress test under frozen rules

(A) Single-country versus multicountry distribution among 45 Correa-labelled studies. (B) Pilot cohort patient counts or selected-isolate counts, according to the frozen metadata status. (C) Locus-screen attrition. (D) Strict endpoint status. Both within-cohort disease contrasts are NOT_COMPUTABLE under the pre-specified joint-callability floor; displayed permutation values are diagnostic only.


## References

1. Tomb JF, White O, Kerlavage AR, Clayton RA, Sutton GG, et al. The complete genome sequence of the gastric pathogen *Helicobacter pylori*. Nature. 1997;388:539–547. doi:10.1038/41483.
2. Saunders NJ, Peden JF, Hood DW, Moxon ER. Simple sequence repeats in the *Helicobacter pylori* genome. Mol Microbiol. 1998;27:1091–1098. doi:10.1046/j.1365-2958.1998.00768.x.
3. Appelmelk BJ, Martin SL, Monteiro MA, Clayton CA, McColm AA, et al. Phase variation in *Helicobacter pylori* lipopolysaccharide due to changes in the lengths of poly(C) tracts in alpha3-fucosyltransferase genes. Infect Immun. 1999;67:5361–5366. doi:10.1128/IAI.67.10.5361-5366.1999.
4. Salaün L, Ayraud S, Saunders NJ. Phase variation mediated niche adaptation during prolonged experimental murine infection with *Helicobacter pylori*. Microbiology. 2005;151:917–923. doi:10.1099/mic.0.27379-0.
5. Moxon R, Bayliss C, Hood D. Bacterial contingency loci: the role of simple sequence DNA repeats in bacterial adaptation. Annu Rev Genet. 2006;40:307–333. doi:10.1146/annurev.genet.40.110405.090442.
6. Gymrek M. A genomic view of short tandem repeats. Curr Opin Genet Dev. 2017;44:9–16. doi:10.1016/j.gde.2017.01.012.
7. Willems T, Zielinski D, Yuan J, Gordon A, Gymrek M, Erlich Y. Genome-wide profiling of heritable and de novo STR variations. Nat Methods. 2017;14:590–592. doi:10.1038/nmeth.4267.
8. Achtman M, Azuma T, Berg DE, Ito Y, Morelli G, et al. Recombination and clonal groupings within *Helicobacter pylori* from different geographical regions. Mol Microbiol. 1999;32:459–470. doi:10.1046/j.1365-2958.1999.01382.x.
9. Falush D, Wirth T, Linz B, Pritchard JK, Stephens M, et al. Traces of human migrations in *Helicobacter pylori* populations. Science. 2003;299:1582–1585. doi:10.1126/science.1080857.
10. Linz B, Balloux F, Moodley Y, Manica A, Liu H, et al. An African origin for the intimate association between humans and *Helicobacter pylori*. Nature. 2007;445:915–918. doi:10.1038/nature05562.
11. Thorell K, Muñoz-Ramírez ZY, Wang D, Sandoval-Motta S, Boscolo Agostini R, et al.; HpGP Research Network. The *Helicobacter pylori* Genome Project: insights into *H. pylori* population structure from analysis of a worldwide collection of complete genomes. Nat Commun. 2023;14:8184. doi:10.1038/s41467-023-43562-y.
12. Correa P. Human gastric carcinogenesis: a multistep and multifactorial process. Cancer Res. 1992;52:6735–6740.
13. Wenger AM, Peluso P, Rowell WJ, Chang PC, Hall RJ, et al. Accurate circular consensus long-read sequencing improves variant detection and assembly of a human genome. Nat Biotechnol. 2019;37:1155–1162. doi:10.1038/s41587-019-0217-9.
14. Ondov BD, Treangen TJ, Melsted P, Mallonee AB, Bergman NH, Koren S, Phillippy AM. Mash: fast genome and metagenome distance estimation using MinHash. Genome Biol. 2016;17:132. doi:10.1186/s13059-016-0997-x.
15. Pazos A, Kodaman N, Piazuelo MB, Romero-Gallo J, Sobota RS, et al. Draft genome sequences of 13 Colombian *Helicobacter pylori* strains isolated from Pacific Coast and Andean residents. Genome Announc. 2017;5:e00113-17. doi:10.1128/genomeA.00113-17.
16. Jackson LK, Potter B, Schneider S, Fitzgibbon M, Blair K, et al. *Helicobacter pylori* diversification during chronic infection within a single host generates sub-populations with distinct phenotypes. PLoS Pathog. 2020;16:e1008686. doi:10.1371/journal.ppat.1008686.
17. Hu L, Zeng X, Ai Q, Liu C, Zhang X, Chen Y, Liu L, Li GQ. Long-read- and short-read-based whole-genome sequencing reveals the antibiotic resistance pattern of *Helicobacter pylori*. Microbiol Spectr. 2023;11:e04522-22. doi:10.1128/spectrum.04522-22.
18. Kong PF, Yan YH, Duan YT, Fang YT, Dou Y, Xu YH, Xu DZ. Comparative genomic analysis of *Helicobacter pylori* isolates from gastric cancer and gastritis in China. BMC Cancer. 2025;25:628. doi:10.1186/s12885-025-13493-6.
19. Zhang X, Liu H, Xu S, Zhang S, Yang T, Lei Z, et al. Within-host diversity and phased variant analysis reveal structures and recombination of *Helicobacter pylori* subpopulations in stomach. GigaScience. 2026;15:giag046. doi:10.1093/gigascience/giag046.


## References

1. Saunders NJ, Peden JF, Hood DW, Moxon ER. Simple sequence repeats in the *Helicobacter pylori* genome. *Mol Microbiol*. 1998;27:1091–1098. PMID: 9570395.

2. Appelmelk BJ, Martin SL, Monteiro MA, et al. Phase variation in *Helicobacter pylori* lipopolysaccharide due to changes in the lengths of poly(C) tracts in alpha3-fucosyltransferase genes. *Infect Immun*. 1999;67:5361–5366. doi:10.1128/IAI.67.10.5361-5366.1999.

3. de Vries N, Duinsbergen D, Kuipers EJ, et al. Transcriptional phase variation of a type III restriction-modification system in *Helicobacter pylori*. *J Bacteriol*. 2002;184:6615–6623. PMCID: PMC135423.

4. Salaün L, Ayraud S, Saunders NJ. Phase variation mediated niche adaptation during prolonged experimental murine infection with *Helicobacter pylori*. *Microbiology*. 2005;151:917–923. doi:10.1099/mic.0.27379-0.

5. Salaün L, Linz B, Suerbaum S, Saunders NJ. The diversity within an expanded and redefined repertoire of phase-variable genes in *Helicobacter pylori*. *Microbiology*. 2004;150:817–830. doi:10.1099/mic.0.26993-0.

6. Falush D, Wirth T, Linz B, et al. Traces of human migrations in *Helicobacter pylori* populations. *Science*. 2003;299:1582–1585. doi:10.1126/science.1080857.

7. Linz B, Balloux F, Moodley Y, et al. An African origin for the intimate association between humans and *Helicobacter pylori*. *Nature*. 2007;445:915–918. doi:10.1038/nature05562.

8. Suerbaum S, Achtman M. *Helicobacter pylori*: recombination, population structure and human migrations. *Int J Med Microbiol*. 2004;294:133–139. doi:10.1016/j.ijmm.2004.06.014.

9. Thorell K, Muñoz-Ramírez ZY, Wang D, et al. The *Helicobacter pylori* Genome Project: insights into *H. pylori* population structure from analysis of a worldwide collection of complete genomes. *Nat Commun*. 2023;14:8184. doi:10.1038/s41467-023-43562-y.

10. Correa P, Piazuelo MB. The gastric precancerous cascade. *J Dig Dis*. 2012;13:2–9. doi:10.1111/j.1751-2980.2011.00550.x.

11. Correa P. *Helicobacter pylori* infection and gastric adenocarcinoma. 2011. PMID: 21857882.

12. Malfertheiner P, Camargo MC, El-Omar E, et al. *Helicobacter pylori* infection. *Nat Rev Dis Primers*. 2023;9:19. doi:10.1038/s41572-023-00431-8.

13. Gymrek M. A genomic view of short tandem repeats. *Curr Opin Genet Dev*. 2017;44:9–16. doi:10.1016/j.gde.2017.01.012.

14. Mousavi N, Shleizer-Burko S, Yanicky R, Gymrek M. Profiling the genome-wide landscape of tandem repeat expansions. *Nucleic Acids Res*. 2019;47:e90. doi:10.1093/nar/gkz501.

15. Microbiology Society. Microbial Genomics: aims and scope; prepare-an-article guidance; author submission checklist. Accessed 24 September 2026.
