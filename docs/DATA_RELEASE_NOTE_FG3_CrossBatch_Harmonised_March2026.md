# FinnGen 3 Cross-Batch Harmonised Proteomics Data Release Note
## Batch 01 + Batch 02 Aggregate — Production Run, February 2026

- **Release Date**: March 2026 (v1.0)
- **Platform**: Olink Explore HT (5K)
- **Batches**: FG3 Batch 01 and FG3 Batch 02 (harmonised aggregate)
- **Input (post-QC)**: Batch 01: 1,957 samples × 5,440 proteins | Batch 02: 2,458 samples × 5,416 proteins
- **Total Measurements (final aggregate)**: 22,303,288 (4,118 samples × 5,416 proteins)
- **Author**: Reza Jabal, PhD (rjabal@broadinstitute.org)
- **Reviewer**: Mitja Kurki, PhD (mkurki@broadinstitute.org)
- **Pipeline**: [https://github.com/FINNGEN/FinnGen3-proteomics-pipeline](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline) v1.7.1

> **Scope**: This release note covers [Steps 06–11](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline#pipeline-steps) of the [FinnGen 3 Olink pipeline](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline) — within-batch normalisation, cross-batch bridge harmonisation, covariate adjustment, phenotype preparation, kinship filtering, and inverse rank normalisation — applied to the post-QC outputs of Batch 01 and Batch 02. Per-batch QC ([Steps 00–05d](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline#pipeline-steps)) is documented in the respective per-batch release notes.

---

## File Descriptions

All files are provided in three formats (RDS, Parquet, TSV) unless otherwise noted. Row identifiers are SampleIDs for per-batch [Step 07](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline#pipeline-steps) matrices and FINNGENIDs for all [Step 09](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline#pipeline-steps)+ matrices and aggregate outputs. All matrices contain 5,416 biological proteins (columns).

### 1. Per-Batch Bridge-Normalised NPX ([Step 07](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline#pipeline-steps) — not covariate-adjusted, not kinship-filtered)

**Files**:
- `npx_bridge_normalised_1957_samples_fg3_batch_01.{rds,parquet,tsv}`
- `npx_bridge_normalised_2458_samples_fg3_batch_02.{rds,parquet,tsv}`

**Description**: Per-batch NPX matrices after within-batch median normalisation ([Step 06](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline#pipeline-steps)) and cross-batch Olink standard bridge normalisation ([Step 07](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline#pipeline-steps)). These are **not covariate-adjusted** and **not kinship-filtered**. Row IDs are the original Olink SampleIDs (not FINNGENIDs). To map SampleIDs to FINNGENIDs, use the per-batch sample mapping file produced in [Step 00](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline#pipeline-steps) (`00_sample_mapping_fg3_batch_01.rds` / `00_sample_mapping_fg3_batch_02.rds`), which contains a `SampleID` → `FINNGENID` lookup for all FinnGen participants. Batch 02 is the reference batch (unchanged during bridge normalisation); Batch 01 has been additively adjusted to align with Batch 02. Use these for analyses requiring the harmonised NPX values before any covariate adjustment or sample removal.

| File | Samples | Proteins | Row identifier |
|---|---|---|---|
| `…_1957_samples_fg3_batch_01` | 1,957 | 5,416 | SampleID |
| `…_2458_samples_fg3_batch_02` | 2,458 | 5,416 | SampleID |

### 2. Per-Batch FINNGENID Pre-Kinship NPX ([Step 09](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline#pipeline-steps) — not kinship-filtered)

**Files**:
- `npx_finngenid_prekinship_1955_samples_fg3_batch_01.{rds,parquet,tsv}`
- `npx_finngenid_prekinship_2448_samples_fg3_batch_02.{rds,parquet,tsv}`

**Description**: Per-batch NPX matrices after covariate adjustment ([Step 08](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline#pipeline-steps)), outlier reconciliation, and FINNGENID conversion ([Step 09](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline#pipeline-steps)). These are **not kinship-filtered** and **not deduplicated** — where a FINNGENID has multiple SampleIDs (technical replicates), all measurements are preserved as separate rows. This is by design: duplicate resolution (selecting the best-quality measurement per FINNGENID) and kinship filtering occur in [Step 10](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline#pipeline-steps). Row IDs are FINNGENIDs (some FINNGENIDs may appear more than once). For Batch 02, the 50 bridging samples have been re-indexed from their EA5_OLI_ SampleIDs to FINNGENIDs using `FG_cross_batch_bridging_metadata.tsv`.

| File | Rows | Unique FINNGENIDs | Duplicate rows | Proteins |
|---|---|---|---|---|
| `…_1955_samples_fg3_batch_01` | 1,955 | 1,933 | 22 | 5,416 |
| `…_2448_samples_fg3_batch_02` | 2,448 | 2,318 | 130 | 5,416 |

### 3. Aggregate Pre-Kinship NPX — Not Batch-Corrected ([Step 09](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline#pipeline-steps))

**Files**: `npx_aggregate_prekinship_4360_samples.{rds,parquet,tsv}`

**Description**: Combined aggregate NPX matrix (Batch 01 + Batch 02) after cross-batch overlap resolution and bridging sample FINNGENID recovery. This matrix is covariate-adjusted (age + sex, [Step 08](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline#pipeline-steps)) and bridge-normalised ([Step 07](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline#pipeline-steps)), but has **not** had the additional per-protein residual batch adjustment applied (`lm(NPX ~ batch)` from [Step 09](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline#pipeline-steps)). Duplicate FINNGENIDs (technical replicates) are preserved. **Not kinship-filtered** and **not rank-normalised**. Row IDs are FINNGENIDs. Use this for analyses where batch-level signal should be preserved or where analysts wish to apply their own batch-correction strategy.

**Dimensions**: **4,360 rows** × 5,416 proteins.

### 3b. Aggregate Pre-Kinship NPX — Batch-Corrected ([Step 09](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline#pipeline-steps))

**Files**: `npx_aggregate_batch_corrected_prekinship_4360_samples.{rds,parquet,tsv}`

**Description**: Same as Section 3 above, but with the additional per-protein residual batch adjustment applied (`lm(NPX ~ batch)`, see Section "Residual Batch Adjustment" below). Duplicate FINNGENIDs (technical replicates) are preserved. **Not kinship-filtered** and **not rank-normalised**. Row IDs are FINNGENIDs.

**Dimensions**: **4,360 rows** × 5,416 proteins.

### 4. Final Aggregate — Kinship-Filtered, Rank-Normalised (pQTL-ready)

**Files**:
- `npx_aggregate_rank_normalised_unrelated_4118_samples.{rds,parquet,tsv}`
- `npx_aggregate_rank_normalised_unrelated_4118_samples.pheno` (PLINK format)

**Description**: The primary analysis-ready dataset for pQTL and GWAS analyses. Contains **4,118 unrelated FinnGen samples** with inverse rank-normalised protein expression values after cross-batch bridge harmonisation, batch correction, kinship filtering, and column-wise inverse rank normalisation. Row IDs are FINNGENIDs. The `.pheno` file is in PLINK phenotype format (FID IID followed by 5,416 protein columns).

**Dimensions**: **4,118 samples** × 5,416 proteins.

### 5. LOD-Filtered High-Quality Assay Subsets

> **Analyst recommendation**: Although this release provides cross-batch normalised data across all 5,416 biological proteins, we recommend that analysts restrict quantitative analyses (e.g. differential expression, association testing) to the **3,075 high-quality assays** where at least 10% of biological samples have NPX values above the per-protein limit of detection (LOD). Proteins that fail this criterion (>90% of samples below LOD) have the vast majority of measurements at or near background noise, which can inflate false-positive rates and reduce statistical power. The full 5,416-protein matrices are provided for completeness and for analysts who wish to apply their own detection thresholds (e.g. a stricter 50% LOD threshold yields 2,102 assays).

**LOD definition**: Per-protein LOD is computed from Negative Control (NC) wells in the raw Olink Explore HT data (Batch 02) as:

`LOD_protein = mean(NC_NPX) + 3 × SD(NC_NPX)`

pooled across all 30 plates (60 NC observations per protein). A protein fails the recommended LOD filter if more than 90% of biological samples have raw NPX < LOD.

**LOD statistics** (5,416 proteins):
- Proteins above LOD (≤90% below): **3,075** (56.8%)
- Proteins below LOD (>90% below): **2,341** (43.2%)
- Proteins with 0% below LOD: 620 (11.4%)
- Median per-protein below-LOD rate: 82.3%

| LOD threshold | Proteins passing | Fraction |
|---|---|---|
| ≤50% below LOD (stringent) | 2,102 | 38.8% |
| ≤75% below LOD | 2,495 | 46.1% |
| ≤90% below LOD (recommended) | **3,075** | **56.8%** |
| ≤95% below LOD (lenient) | 3,564 | 65.8% |

**Files**:

| File | Description |
|---|---|
| `assay_list_above_lod_3075.txt` | Plain-text list of 3,075 protein names that pass the 90% LOD filter (one per line) |
| `lod_nc_summary_batch02.tsv` | Full per-protein LOD summary: Assay, OlinkID, n_bio, n_below_lod, pct_below_lod, LOD, nc_mean, nc_sd. Analysts can use this table to apply any custom threshold |
| `npx_aggregate_above_lod_prekinship_4360_samples_3075_proteins.{rds,parquet,tsv}` | Pre-kinship aggregate (not batch-corrected) sliced to 3,075 high-quality proteins (4,360 × 3,075). Not kinship-filtered, not rank-normalised. Row IDs: FINNGENIDs |
| `npx_bridge_normalised_above_lod_1957_samples_3075_proteins_fg3_batch_01.{rds,parquet,tsv}` | Batch 01 bridge-normalised NPX, sliced to 3,075 high-quality proteins. Not covariate-adjusted, not kinship-filtered. Row IDs: SampleIDs |
| `npx_bridge_normalised_above_lod_2458_samples_3075_proteins_fg3_batch_02.{rds,parquet,tsv}` | Batch 02 bridge-normalised NPX, sliced to 3,075 high-quality proteins. Not covariate-adjusted, not kinship-filtered. Row IDs: SampleIDs |

### 6. Auxiliary Files

| File | Description |
|---|---|
| `sample_list_aggregate_unrelated_4118.txt` | Plain-text list of 4,118 FINNGENIDs in the final unrelated aggregate (one per line) |
| `protein_list_5416.txt` | Plain-text list of 5,416 protein identifiers (one per line), matching column order of all matrices |
| `batch_r2_diagnostics.tsv` | Per-protein batch R² from `lm(NPX ~ batch)` on the pre-correction combined matrix |
| `comprehensive_outliers_list_fg3_batch_01.tsv` | QC outliers for Batch 01 ([Step 05d](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline#pipeline-steps)) with method flags and metrics |
| `comprehensive_outliers_list_fg3_batch_02.tsv` | QC outliers for Batch 02 ([Step 05d](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline#pipeline-steps)) with method flags and metrics |
| `removed_221_breakdown_summary.tsv` | Tabular breakdown of the 221 samples removed between QC-passed and final aggregate |
| `removed_221_full_id_details.tsv` | Full SampleID/FINNGENID-level detail for all 221 removed samples |

---

## Input Summary — QC-Passed Starting Points

This processing pipeline receives as input the per-batch post-QC NPX matrices produced by [Steps 00–05d](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline#pipeline-steps).

| Batch | Samples | Proteins | Source |
|---|---|---|---|
| Batch 01 | 1,957 | 5,440 | [Step 05d](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline#pipeline-steps) QC-passed matrix |
| Batch 02 | 2,458 | 5,416 | [Step 05d](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline#pipeline-steps) QC-passed matrix |

**Protein panel**: Batch 01 contains 5,440 proteins (including 24 control probes); Batch 02 contains 5,416 biological proteins. All cross-batch and aggregate outputs use the **intersection of 5,416 common proteins**, excluding the 24 control probes present in Batch 01.

### Sex Outlier Detection Mode

The production run used `sex_outlier_mode: strict_only` for both batches. Under this mode:

- **TIER 1 (strict mismatches)**: Samples where `predicted_sex ≠ genetic_sex` are flagged **and removed** from the QC-passed dataset.
- **TIER 2 (threshold-based outliers)**: Samples with borderline sex prediction probabilities that still match genetic sex are flagged but **kept** in the dataset.

| Batch | Youden's J threshold | Strict mismatches (removed) | Threshold-based outliers (kept) | Total flagged |
|---|---|---|---|---|
| Batch 01 | 0.50 (Platt-scaled) | 45 | 0 | 45 |
| Batch 02 | 0.63 | 17 | 6 | 23 |

This differs from the December 2025 Batch 02 QC release (v.2.1) which used `sex_outlier_mode: all`, where both strict mismatches and threshold-based outliers would be removed. The `strict_only` mode was chosen for the cross-batch harmonisation to maximise sample retention whilst removing only high-confidence sex-provenance errors.

---

## Within-Batch Normalisation ([Step 06](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline#pipeline-steps))

**Purpose**: Equalise the median protein expression level across all samples within each batch to remove systematic technical shifts in overall signal intensity. Performed independently per batch.

**Method**: Median normalisation — for each protein, the sample-level median NPX is subtracted and the global batch median added back.

**Results**:

| Metric | Batch 01 (before → after) | Batch 02 (before → after) |
|---|---|---|
| Mean SD | 1.151 → 0.293 | 1.128 → 0.263 |
| Mean MAD | 0.876 → 0.220 | 0.896 → 0.201 |
| Mean IQR | 1.271 → 0.320 | 1.262 → 0.287 |
| Mean SD reduction | **74.5%** | **76.7%** |

No samples are removed at this step.

---

## Cross-Batch Bridge Harmonisation ([Step 07](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline#pipeline-steps))

### Bridge Sample Design

| Parameter | Value |
|---|---|
| Bridge FINNGENIDs identified (labelled as `Bridging` in both batches) | 24 |
| Bridge samples in Batch 01 QC-passed matrix | 23 (1 removed: SampleID `52291062464138167509` failed Batch 01 QC) |
| Bridge samples in Batch 02 QC-passed matrix | 23 (1 removed: failed Batch 02 QC) |
| **Matched bridge pairs** (present in both QC-passed matrices) | **22** |
| Common proteins | 5,416 |
| **Reference batch** | **Batch 02** (unchanged) |

**Note**: 30 FINNGENIDs are recorded as being in both batches in `FG_cross_batch_bridging_metadata.tsv`, but only 24 of these were labelled as `sample_type = "Bridging"` in the Batch 01 sample mapping. The remaining 6 were classified as regular FinnGen participants in Batch 01 and were not used for bridge normalisation. Of the 24 labelled bridging samples, 2 different individuals failed per-batch QC (one in each batch), leaving 22 matched pairs with NPX data in both QC-passed matrices.

### Harmonisation Method — Olink Standard Bridge Normalisation

For each protein `p` and each of the 22 matched bridge sample pairs `k`, the pairwise difference is computed:

`diff_k,p = NPX_reference[k, p] − NPX_target[k, p]`

The per-protein adjustment factor is the median across all 22 pairs:

`Adj_factor[p] = median(diff_1,p, diff_2,p, …, diff_22,p)`

The target batch (Batch 01) is then additively corrected per protein: `NPX_Batch01_harmonised[p] = NPX_Batch01[p] + Adj_factor[p]`. Batch 02 (reference) is unchanged.

### Method Comparison

Three cross-batch harmonisation methods were evaluated on the same input data. **ComBat** (empirical Bayes batch correction from the `sva` R package) estimates and removes batch effects using a parametric or non-parametric empirical Bayes framework; it was included as a comparison method but adjusts both batches. **Cross-batch median** applies a global median alignment between batches.

| Method | Batch 01 CV reduction | Batch 02 (reference) |
|---|---|---|
| **Bridge (Olink standard)** | **42.2%** | 0% (unchanged) |
| ComBat (empirical Bayes) | 14.8% | −21.7% (worsened) |
| Cross-batch median | −518% (worsened) | 0% (unchanged) |

Bridge normalisation was selected because it achieves the largest CV reduction in the adjusted batch while leaving the reference batch completely unchanged. No samples are removed at this step.

---

## Harmonisation Key Performance Indicator (KPI) Evaluation ([Step 07b](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline#pipeline-steps))

A post-hoc quantitative evaluation of harmonisation quality was performed using the 22 matched bridge pairs. The KPI framework compares the pre-harmonisation ([Step 05d](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline#pipeline-steps)) and post-harmonisation ([Step 07](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline#pipeline-steps)) matrices across multiple distance, concordance, and biological preservation metrics to objectively assess whether bridge normalisation reduced inter-batch technical variation without degrading biological signal.

### Bridge Pair Convergence — Passing Metrics

| Metric | Pre | Post | Result |
|---|---|---|---|
| Median bridge Euclidean distance | 50.13 | 38.84 | **PASS** (decreased) |
| Median bridge Mahalanobis distance | 11.14 | 7.16 | **PASS** (decreased) |
| Median bridge MAD-whitened distance (primary) | 8.73 | 6.46 | **PASS** (26% reduction) |
| Wilcoxon signed-rank p-value | — | 3 × 10⁻⁶ | **PASS** (< 0.05) |
### Noise Floor Analysis

The noise floor establishes a lower bound on the distance we can expect between two measurements of the same individual, even without any inter-batch technical variation. It is computed from within-batch replicate pairs — individuals who were measured more than once within the same batch — using the median absolute deviation (MAD) of their pairwise Euclidean distances in PC space. If bridge normalisation is working correctly, the post-harmonisation bridge pair distance should approach this within-batch replicate noise, since a perfectly harmonised bridge pair should look no more different than two replicate measurements of the same sample.

| Metric | Value |
|---|---|
| Within-batch replicate noise (MAD) | 7.24 |
| Post-harmonisation bridge pair distance (MAD-whitened) | 6.46 |
| **Noise floor ratio** (bridge post / replicate noise) | **0.89** |

A noise floor ratio close to 1.0 indicates that post-harmonisation bridge pair distances are approximately equal to the irreducible technical measurement noise. The observed ratio of 0.89 confirms that the bridge normalisation has brought cross-batch differences down to (or slightly below) the level of within-batch replicate variation — the best achievable result.

### Biological Signal Preservation

A critical requirement of any batch harmonisation method is that it removes technical inter-batch variation without degrading genuine biological signal. To verify this, two biological metrics were evaluated before and after harmonisation:

1. **Sex variance fraction** — the median fraction of per-protein expression variance attributable to sex (from a per-protein ANOVA decomposing variance into batch, sex, and residual components). If harmonisation were inadvertently removing sex-related biological signal, this fraction would decrease.

2. **Sex prediction accuracy** — the accuracy of predicting genetic sex from the protein expression matrix using a simple classifier. A drop in accuracy after harmonisation would indicate that biologically meaningful signal had been distorted.

| Metric | Pre-harmonisation | Post-harmonisation | Assessment |
|---|---|---|---|
| Sex variance fraction (median across 5,416 proteins) | 6 × 10⁻⁴ | 6 × 10⁻⁴ | **Stable** — no biological signal loss |
| Sex prediction accuracy | 0.596 | 0.596 | **Stable** — classification unaffected |

Both metrics are unchanged, confirming that bridge normalisation preserved biologically relevant protein expression patterns while removing inter-batch technical variation.

---

## Covariate Adjustment ([Step 08](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline#pipeline-steps))

Covariate adjustment was performed independently per batch using linear regression (`NPX ~ age + sex`) on samples with complete covariate data. Model residuals (plus protein mean) replace the original NPX values.

**Note**: Proteomic PCs (pPC1–10) were evaluated and visualised but **NOT** used in the adjustment.

| Parameter | Batch 01 | Batch 02 |
|---|---|---|
| Input samples | 1,957 | 2,458 |
| Samples with complete covariates | 1,447 (73.9%) | 1,545 (62.9%) |
| Output samples | 1,957 | 2,458 |


---

## Phenotype Preparation and Sample Reconciliation ([Step 09](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline#pipeline-steps))

### Outlier Reconciliation and FINNGENID Conversion

| Batch | Input | Outlier reconciliation | After conversion |
|---|---|---|---|
| Batch 01 | 1,957 | −2 | 1,955 |
| Batch 02 | 2,458 | −10 | 2,448 |

**Bridging sample FINNGENID recovery**: The Batch 02 QC-passed matrix contains 50 cross-batch bridging samples measured under EA5 tube identifiers (EA5_OLI_ SampleIDs). The [FinnGen 3 Olink pipeline](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline)'s [Step 00](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline#pipeline-steps) sample mapping did not resolve these EA5-format SampleIDs to FINNGENIDs automatically. Because these are genuine FinnGen participants with known FINNGENIDs (recorded in `FG_cross_batch_bridging_metadata.tsv`), we recovered their identities by mapping each EA5_OLI_ SampleID against the `PSEUDO_ID_Batch_00` column in that file and re-indexed them by FINNGENID. All 50 bridging samples are now present in the FINNGENID-indexed pre-kinship matrices under their correct FINNGENIDs.

### Duplicate FINNGENID Selection

Some FINNGENIDs are associated with multiple Olink SampleIDs (technical replicates — the same individual measured more than once). To produce a matrix with exactly one measurement per individual, the pipeline selects the single highest-quality SampleID for each duplicated FINNGENID using a deterministic hierarchical criterion: (1) lowest missing rate, (2) highest number of detected proteins, (3) lowest SD NPX, (4) lexicographically first SampleID. The remaining lower-quality replicates are excluded.

| Batch | Duplicate groups | Samples removed |
|---|---|---|
| Batch 01 | 22 | 22 |
| Batch 02 | 130 | 130 |

### Cross-Batch Overlap Resolution

When the two per-batch FINNGENID matrices are combined into a single aggregate, some FINNGENIDs appear in both Batch 01 and Batch 02. These overlapping individuals fall into two distinct categories:

1. **FinnGen participants measured independently in both batches** (35 individuals) — these are regular FinnGen samples whose plasma was processed in both experimental runs. They are not bridging samples; their presence in both batches is coincidental (e.g. the same biobank participant's sample was included in two separate Olink runs).

2. **Bridging samples measured in both batches by design** (6 individuals) — these are the subset of the 30 cross-batch bridging FINNGENIDs (from `FG_cross_batch_bridging_metadata.tsv`) whose Batch 01 SampleIDs were classified as regular FinnGen samples (`sample_type = "FinnGen"`) rather than `"Bridging"` in the Batch 01 sample mapping. They were not used for bridge normalisation in [Step 07](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline#pipeline-steps), but they do have measurements in both batches.

For all 41 overlapping FINNGENIDs, the **Batch 02 measurement is retained** in the aggregate (Batch 02 is the reference batch for bridge normalisation); the corresponding Batch 01 measurement is excluded. This prevents duplicate rows for the same individual in the aggregate matrix.

### Residual Batch Adjustment

The primary inter-batch harmonisation is performed by bridge normalisation in [Step 07](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline#pipeline-steps). However, after aggregating the two batches, a small residual linear batch signal may remain for some proteins. To address this, an additional per-protein linear adjustment was applied on the combined aggregate matrix:

For each protein: `lm(NPX ~ batch)` — the model residuals (plus the protein grand mean) replace the original values. This removes only the linear component of any remaining batch signal; it does not replace or repeat the bridge normalisation, which operates at the individual bridge-sample level.

**Batch R² statistics** (per-protein R² from the batch-effect model before adjustment):

| Statistic | Value |
|---|---|
| Median batch R² | 0.020 |
| Mean batch R² | 0.057 |
| Proteins with R² > 0.10 | 18.5% (1,002/5,416) |
| Proteins with R² > 0.50 | 0.7% (38/5,416) |

The low median R² (0.02) confirms that bridge normalisation ([Step 07](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline#pipeline-steps)) already removed the majority of the inter-batch technical signal; this additional adjustment addresses the residual tail.

---

## Kinship Filtering ([Step 10](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline#pipeline-steps))

**Threshold**: KING kinship coefficient ≥ 0.0884 (2nd-degree relationships and closer: parent-offspring, full siblings, half-siblings, grandparent-grandchild).

**Within-batch kinship filtering** ([Step 10](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline#pipeline-steps) of the pipeline):

| Batch | Related pairs | Samples removed | Unrelated retained |
|---|---|---|---|
| Batch 01 | 16 (13 first-degree + 3 second-degree) | 16 | 1,882 |
| Batch 02 | 7 (4 first-degree + 3 second-degree) | 7 | 2,262 |

**Cross-batch kinship filtering** (post-hoc, applied to the combined aggregate):

After combining the per-batch kinship-filtered matrices, 27 cross-batch related pairs (Kinship ≥ 0.0884) were identified — these are pairs where one member is in Batch 01 and the other in Batch 02, which were not assessed during per-batch filtering (9 parent-offspring, 5 full-sibling, 13 second-degree). For each pair, the Batch 01 sample was removed (retaining the Batch 02 measurement, consistent with the reference batch convention). Of the 27 pairs, one Batch 01 individual appeared in two separate related pairs; removing that single sample resolved both pairs simultaneously, yielding 26 unique removals from 27 pairs. All 26 removed samples originated from Batch 01.

| Stage | Samples removed | Unrelated retained |
|---|---|---|
| Within-batch (B01 + B02) | 23 | 4,144 |
| Cross-batch | 26 | **4,118** |
| **Total kinship removals** | **49** | **4,118** |

---

## Rank Normalisation ([Step 11](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline#pipeline-steps))

Inverse rank normalisation (IRN) applied column-wise per protein. Each protein's values are replaced by the standard normal quantile of their within-protein rank.

| Metric | Batch 01 | Batch 02 |
|---|---|---|
| Samples | 1,882 | 2,262 |
| Mean skewness (before → after) | 0.944 → 0.000 | 0.815 → 0.000 |
| Mean kurtosis (before → after) | 9.571 → 2.984 | 7.196 → 2.986 |

**Aggregate**: 4,118 samples × 5,416 proteins = 22,303,288 total measurements.

---

## Disease Group Summary

Summary of disease groups across the two key release matrices. Groups are assigned based on available metadata; samples without disease information but confirmed as FinnGen bridging samples are labelled "Bridging".

### 1. Pre-Kinship Aggregate (4,360 rows, duplicates preserved)

This matrix contains all measurements before kinship filtering and duplicate resolution. Counts represent rows (SampleIDs), not unique individuals.

| Disease Group | Batch 01 | Batch 02 | Aggregate |
|---|---|---|---|
| Blood_Donor | 1,108 | 394 | 1,480 |
| Kidney | 234 | 657 | 886 |
| Variant_carriers | 380 | 11 | 383 |
| Parkinsons | 0 | 335 | 335 |
| Rheuma | 0 | 260 | 259 |
| IBD | 228 | 8 | 229 |
| AMD | 0 | 162 | 162 |
| Chromosomal_Abnormalities | 0 | 153 | 153 |
| Metabolic | 0 | 133 | 133 |
| MFGE8 | 0 | 117 | 117 |
| F64 | 0 | 61 | 61 |
| Kids | 0 | 53 | 53 |
| Pulmo | 0 | 46 | 46 |
| Bridging | 0 | 43 | 43 |
| Multiple_Diseases | 5 | 15 | 20 |
| **Total** | **1,955** | **2,448** | **4,360** |

**Note**: The Aggregate column differs from Batch 01 + Batch 02 because cross-batch overlap resolution removes 41 Batch 01 rows (FINNGENIDs present in both batches keep the Batch 02 measurement). The "Bridging" group consists of the 43 cross-batch bridge samples in Batch 02 that do not have FinnGen phenotype metadata.

### 2. Final pQTL Aggregate (4,118 unique FINNGENIDs)

This matrix is deduplicated (one measurement per FINNGENID) and kinship-filtered (unrelated individuals only). Batch origin is determined by cross-batch overlap resolution: FINNGENIDs present in both batches use the Batch 02 measurement.

| Disease Group | Batch 01 | Batch 02 | Total | % |
|---|---|---|---|---|
| Blood_Donor | 1,057 | 391 | 1,448 | 35.2% |
| Kidney | 228 | 653 | 881 | 21.4% |
| Variant_carriers | 363 | 10 | 373 | 9.1% |
| Parkinsons | 0 | 334 | 334 | 8.1% |
| Rheuma | 0 | 235 | 235 | 5.7% |
| IBD | 197 | 7 | 204 | 5.0% |
| AMD | 0 | 160 | 160 | 3.9% |
| Chromosomal_Abnormalities | 0 | 149 | 149 | 3.6% |
| MFGE8 | 0 | 116 | 116 | 2.8% |
| Metabolic | 0 | 67 | 67 | 1.6% |
| F64 | 0 | 60 | 60 | 1.5% |
| Kids | 0 | 53 | 53 | 1.3% |
| Pulmo | 0 | 23 | 23 | 0.6% |
| Multiple_Diseases | 5 | 10 | 15 | 0.4% |
| **Total** | **1,850** | **2,268** | **4,118** | **100%** |

---

## Overall Sample Flow Summary

```
═══════════════════════════════════════════════════════════════════════
 BATCH 01 post-QC ([Step 05d](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline#pipeline-steps)):     1,957 samples × 5,440 proteins
 BATCH 02 post-QC ([Step 05d](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline#pipeline-steps)):     2,458 samples × 5,416 proteins
 Combined QC-passed (union):       4,415 samples
═══════════════════════════════════════════════════════════════════════

  ↓ [Step 06](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline#pipeline-steps): Within-batch median normalisation   → no sample change
  ↓ [Step 07](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline#pipeline-steps): Cross-batch bridge harmonisation    → no sample change
  ↓ [Step 08](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline#pipeline-steps): Covariate adjustment (age, sex)     → no sample change

  ↓ STEP 09: Phenotype preparation & reconciliation
     Batch 01 — Outlier reconciliation:               -2
     Batch 01 — Duplicate FINNGENID selection:       -22
     Batch 01 — Cross-batch overlap (→ Batch 02):   -41
     Batch 02 — Outlier reconciliation:              -10
     Batch 02 — Bridging EA5→FINNGENID recovery:     +0 (50 re-indexed, not lost)
     Batch 02 — Duplicate FINNGENID selection:      -130
     Net removed in [Step 09](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline#pipeline-steps):                        -205
     Pre-kinship aggregate (incl. bridging):    4,360 rows (4,210 unique FINNGENIDs)

  ↓ STEP 10: Within-batch kinship filtering (KING ≥ 0.0884)
     Batch 01 — Related samples removed:             -16
     Batch 02 — Related samples removed:              -7
     After within-batch kinship:                 4,144 samples

  ↓ Cross-batch kinship filtering (post-hoc, KING ≥ 0.0884)
     Cross-batch related pairs identified:            27
     Unique Batch 01 samples removed:               -26  (1 individual in 2 pairs)
     After cross-batch kinship:                  4,118 samples

  ↓ [Step 11](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline#pipeline-steps): Inverse rank normalisation          → no sample change

═════════════════════════════════════════════════════════════════════
 AGGREGATE FINAL (batch-corrected, fully unrelated,
 rank-normalised):                4,118 samples × 5,416 proteins
═════════════════════════════════════════════════════════════════════
 Total kinship removals:                      49 (23 within + 26 cross)
 Total samples removed (QCed → final):       297 (6.7% of 4,415)
 Final retention rate:                       93.3% (4,118 / 4,415)
═════════════════════════════════════════════════════════════════════
```

---

## Notes

1. **Protein panel**: The cross-batch intersection is 5,416 biological proteins. The 24 control probes in Batch 01 are excluded from all aggregate outputs.

2. **Sex outlier mode**: The production run used `sex_outlier_mode: strict_only`. Only TIER 1 strict mismatches (predicted_sex ≠ genetic_sex) were removed during per-batch QC. TIER 2 threshold-based outliers were flagged but retained in the dataset. This differs from the December 2025 Batch 02 release (which used mode `all`). The Batch 01 Youden's J threshold was 0.50 (Platt-scaled); Batch 02 was 0.63.

3. **Kinship filtering**: Kinship filtering was applied in two stages using the KING kinship coefficient (threshold ≥0.0884, capturing 2nd-degree relationships and closer: parent-offspring, full siblings, half-siblings, grandparent-grandchild). First, within-batch filtering ([Step 10](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline#pipeline-steps)) removed 23 related pairs (16 from Batch 01, 7 from Batch 02). Second, a post-hoc cross-batch kinship filter was applied to the combined aggregate using `finngen_R13.kin0`, identifying 27 cross-batch related pairs (9 parent-offspring, 5 full-sibling, 13 second-degree) and removing 26 Batch 01 samples (Batch 02 measurements retained). One Batch 01 individual appeared in two separate cross-batch pairs; removing that single individual resolved both pairs simultaneously, hence 27 pairs are resolved by 26 removals. The final 4,118-sample aggregate contains **zero related pairs** at Kinship ≥0.0884.

4. **Duplicate FINNGENID selection**: Where multiple SampleIDs map to the same FINNGENID (technical replicates), only the highest-quality measurement is kept (lowest missing rate → most detected proteins → lowest SD NPX → lexicographic SampleID). Full details in `removed_221_full_id_details.tsv`.

5. **Covariate adjustment coverage**: Samples missing covariate data retain their unadjusted NPX values (no covariate imputation). Complete-covariate fractions: 73.9% (Batch 01) and 62.9% (Batch 02).

6. **Bridge sample recovery**: The 50 bridging samples in Batch 02 (EA5_OLI_ SampleIDs) have been re-indexed to FINNGENIDs using `FG_cross_batch_bridging_metadata.tsv` and are present in the pre-kinship matrices alongside all other measurements. Duplicate FINNGENIDs (technical replicates, including bridging samples that also have P-format SampleIDs) are preserved in the pre-kinship matrices; resolution to one measurement per individual occurs in [Step 10](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline#pipeline-steps) (kinship filtering).

7. **Batch R² interpretation**: Median R² of 0.02 indicates bridge normalisation ([Step 07](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline#pipeline-steps)) already removed the majority of inter-batch signal; the residual batch adjustment in [Step 09](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline#pipeline-steps) addresses the remaining tail. The 0.7% of proteins with R² > 0.50 warrant targeted review.

8. **Pre-kinship aggregate count**: The pre-kinship aggregate matrices contain 4,360 rows (4,210 unique FINNGENIDs, 150 duplicate rows from technical replicates). Batch 01 contributes 1,955 rows (1,933 unique FINNGENIDs) and after cross-batch overlap resolution the remaining B01 contribution plus full B02 (2,448 rows, 2,318 unique FINNGENIDs) yields the 4,360-row aggregate. Duplicate FINNGENID resolution is performed in [Step 10](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline#pipeline-steps).

9. **Missing data**: Batch 01 final rank-normalised matrix has 0.01% missing (1,280 values). Batch 02 has 0% missing. Missing values in PLINK files are represented as `NA`.

10. **Cross-batch overlap**: 41 FINNGENIDs present in both batches (35 FinnGen participants + 6 bridging samples) are resolved to Batch 02 in all aggregate outputs.

11. **Limit-of-detection (LOD) filtering recommendation**: Of the 5,416 proteins in the release, 3,075 (56.8%) have ≤90% of biological samples with raw NPX below the per-protein LOD (computed as mean(NC_NPX) + 3×SD(NC_NPX) from Batch 02 Negative Control wells). The remaining 2,341 proteins have >90% of measurements at or near background and should be treated with caution in quantitative analyses. LOD-filtered subsets of the bridge-normalised and pre-kinship aggregate matrices are provided (see Section 5). The full per-protein LOD summary (`lod_nc_summary_batch02.tsv`) is included so analysts can apply a stricter threshold (e.g. 50% yields 2,102 assays) if required.

---

**Pipeline**: [https://github.com/FINNGEN/FinnGen3-proteomics-pipeline](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline) v1.7.1 (February 2026)
**Production Run**: `output_crossBatch_normalised_pQTL_production_02_26_26` (executed 16–18 February 2026)
**Data Release**: FG3 Batch 01 + Batch 02 Cross-Batch Harmonised Aggregate
**Final Sample Count**: 4,118 samples × 5,416 proteins (batch-corrected, kinship-filtered, rank-normalised)
**Prepared by**: Reza Jabal, PhD — Broad Institute of MIT and Harvard
