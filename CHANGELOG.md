# [1.7.0](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline/compare/v1.6.1...v1.7.0) (2026-02-13)


### Features

* **07b,07:** KPI dashboard enhancements, bridge sample exclusion, and complete-case PCA ([4b96117](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline/commit/4b96117798f8aa3ea0c996e5683f081d308a6f85))

## [1.6.1](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline/compare/v1.6.0...v1.6.1) (2026-02-12)


### Bug Fixes

* fixed isolated bugs in the cross-batch harmonisation KPI evaluation script ([2a9ad03](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline/commit/2a9ad032404dad34f3eae98f3263e9555bf85b88))

# [1.6.0](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline/compare/v1.5.0...v1.6.0) (2026-02-12)


### Features

* **pipeline:** add cross-batch KPI evaluation (step 07b) and enforce QC-passed input for bridge normalisation ([538cb9f](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline/commit/538cb9f46c1421f99883fb43212816c59f23852b))

# [1.5.0](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline/compare/v1.4.0...v1.5.0) (2026-02-05)


### Features

* **sex_outliers:** add configurable sex outlier handling mode ([2ae29f5](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline/commit/2ae29f5365e789558e0c45983d02a2bcc3182df5))

# [1.4.0](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline/compare/v1.3.1...v1.4.0) (2026-01-30)


### Features

* **bridge_normalization:** implement Olink standard bridge normalization with enhanced QC visualizations ([13ebd9c](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline/commit/13ebd9c061ed0292f50cd945dbe56f3c3c972706)), closes [#2166](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline/issues/2166) [#67A9](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline/issues/67A9) [#F4A582](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline/issues/F4A582) [#B2182](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline/issues/B2182)

## [1.3.1](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline/compare/v1.3.0...v1.3.1) (2026-01-29)


### Bug Fixes

* **pipeline:** critical fix for cross-batch normalisation - use additive adjustment ([1f7241f](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline/commit/1f7241f3bce1087d2030b4575849ffaa94b3aebe))

# [1.3.0](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline/compare/v1.2.1...v1.3.0) (2026-01-27)


### Features

* **sex_outliers:** implement conditional Platt scaling for sub-optimal separation ([a3e9819](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline/commit/a3e981937393b883f9916b22e7f0e15f6f045cab))

## [1.2.1](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline/compare/v1.2.0...v1.2.1) (2026-01-26)


### Bug Fixes

* **pipeline:** implement flexible normalization fallback and metadata handling ([8190bca](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline/commit/8190bcaa10cedc7c14a582c01a0609c7298d1b19))

# [1.2.0](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline/compare/v1.1.0...v1.2.0) (2026-01-23)


### Features

* **pipeline:** make covariate adjustment config-driven with age and sex as defaults ([e18afdc](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline/commit/e18afdcb9dad2dba66f66789599324ae445f59c7))

## [1.0.1](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline/compare/v1.0.0...v1.0.1) (2026-01-23)


### Bug Fixes

* **pipeline:** fix step 05c skip logic to prevent execution when disabled ([95c688f](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline/commit/95c688f2cb10d35270ea2fd2479cbdb4c2d95ae2))

# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

* **cross_batch_kpi:** Post-hoc cross-batch harmonisation KPI evaluation (Step 07b)
  * New pipeline step `07b_cross_batch_harmonisation_kpi.R` runs after bridge normalisation
  * Fixed PCA basis (avoids basis leakage) for all distance computations
  * Bridge pair collapse metrics: Euclidean, Mahalanobis, ICC, Lin's CCC
  * Paired log-ratio effect size with 10,000-resample bootstrap 95% CI
  * Within-batch replicate noise floor ratio (23 B1 + 139 B2 replicate pairs)
  * Rank-1 k-NN pair identification rate
  * Batch separability: Silhouette score, kBET acceptance, LISI mixing score
  * Biology preservation: sex prediction accuracy, variance decomposition (batch vs sex)
  * Outlier detection: Mahalanobis flagging, per-protein concordance profiling
  * Multi-panel PDF dashboard with 8+ diagnostic panels
  * KPI summary table with pass/warn/fail classification
  * New config parameter: `parameters.bridge_normalization.run_kpi_evaluation`
  * New R package dependencies: irr, DescTools, FNN, kBET, lisi

* **sex_outliers:** Configurable sex outlier handling mode
  * New configuration parameter `parameters.outliers.sex_outlier_mode` with two options:
    * `"all"` (default): Remove both TIER 1 strict mismatches AND TIER 2 threshold-based outliers
    * `"strict_only"`: Remove only TIER 1 strict mismatches, keep TIER 2 threshold-based outliers
  * Allows analysts to retain borderline cases that may represent valid biological variation
  * Updated scripts: `04_sex_outliers.R`, `05d_qc_comprehensive_report.R`
  * Configuration-aware logging in step 04 and step 05d

* **documentation:** Comprehensive sex outlier detection algorithm specification
  * Created `SEX_OUTLIER_DETECTION_ALGORITHM.md` (untracked development document)
  * Corrected explanation of TIER 2 outliers: they exist in the gap between 0.5 and Youden's J
  * Critical insight: which sex can be TIER 2 outliers depends on Youden's J position
  * When Youden > 0.5: only females can be TIER 2 outliers
  * When Youden < 0.5: only males can be TIER 2 outliers
  * Updated README.md with corrected two-tier QC approach explanation

* **bridge_normalization:** Olink standard bridge normalization implementation
  * Implemented pairwise difference calculation for bridge samples (same FINNGENID across batches)
  * Reference batch remains unchanged (Olink standard approach)
  * Only non-reference batch is adjusted using median of pairwise differences
  * Aligns with official OlinkAnalyze package methodology
  * Added `parameters.bridge_normalization.reference_batch` config option
  * Added `parameters.bridge_normalization.method` config option ("olink_standard" or "combined_reference")
  * Added `parameters.bridge_normalization.bridgeability_proteins` for per-protein QC plots

* **bridge_normalization:** Enhanced PCA visualizations
  * Added silhouette score calculations for clustering quality assessment

* **bridge_normalization:** Reintroduced Coefficient of Variation (CV) metric
  * CV (SD/mean) is now the primary metric for normalization method selection
  * Aligns with Olink standard evaluation approach
  * CV reduction percentage calculated for all normalization methods
  * Best method automatically selected based on highest CV reduction
  * CV metrics logged for both batches before and after normalization

* **bridge_normalization:** Enhanced distribution plots
  * Added kurtosis and skewness measures for before/after normalization
  * Distribution shape metrics displayed in plot subtitles
  * Applied to bridge, median, and ComBat normalization effect plots

* **bridge_normalization:** Improved bridgeability assessment
  * Replaced rudimentary sample counts panel with distribution plots
  * NPX value distributions shown across batches before/after normalization
  * Pairwise t-test p-values included for statistical comparison
  * Per-protein bridgeability plots configurable via YAML
  * Enhanced correlation and effect size metrics

* **README:** Comprehensive documentation updates
  * Detailed explanation of Olink standard bridge normalization algorithm
  * Documentation of pairwise difference calculation approach
  * Reference batch strategy explained
  * Enhanced visualization features documented
  * Evaluation metrics section added (CV, SD, MAD, IQR, kurtosis, skewness, silhouette)
  * Configuration options for bridge normalization detailed

### Changed

* **pipeline:** Make covariate adjustment config-driven with age and sex as defaults
  * Modified step 08 to read covariates from config.yaml as a list
  * Changed default behavior to only adjust for age and sex (removed BMI and smoking from defaults)
  * Added support for `covariates_to_adjust` list in config.yaml
  * Maintained backward compatibility with legacy boolean flags
  * Added validation to ensure only valid covariates (age, sex, bmi, smoking) are used
  * Improved logging to show which covariates are being adjusted

* **config.yaml.template:** Enhanced bridge normalization configuration section
  * Added dedicated `parameters.bridge_normalization` section
  * Reference batch selection (null for auto-selection)
  * Method selection (olink_standard vs combined_reference)
  * Bridgeability proteins list for detailed QC

### Fixed

* **bridge_normalization:** Validation and error handling improvements
  * Added validation for bridge_mapping sample pairing
  * Ensured bridge matrices match expected dimensions
  * Improved error messages for debugging
  * Fixed nested conditional logic for method selection

## [1.0.0] (2026-01-23)

### Bug Fixes

* docker tag error fixed ([6519380](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline/commit/6519380a0d4bd34e386a848fb6f1e4af60c7c03d))
* **pipeline:** make Youden J detection failure-proof and remove hardcoded batch references ([93bf990](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline/commit/93bf99015082a2297e79fb78aa26da20ae33234f))
* **release:** Fix the GitHub Actions npm cache ([c22505f](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline/commit/c22505f50f1a3609d9577e805f8028f73d75b5f7))
* Semantic Release tag detection — enhanced fix ([41c64ff](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline/commit/41c64ff34d748280945fc1a86e90585ecdfe5e1b))
* Semantic Release tag error fixed ([9cbf14b](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline/commit/9cbf14b2002a2dfd171e4d4b735d652c61b5b1db))

### Features

* **release:** implement automated release pipeline with semantic-release ([e84360c](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline/commit/e84360cc11c47a16a880e3d21019af3c705ae76d))
