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
