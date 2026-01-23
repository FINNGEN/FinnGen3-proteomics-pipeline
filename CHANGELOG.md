## [1.0.1](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline/compare/v1.0.0...v1.0.1) (2026-01-23)


### Bug Fixes

* **pipeline:** fix step 05c skip logic to prevent execution when disabled ([95c688f](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline/commit/95c688f2cb10d35270ea2fd2479cbdb4c2d95ae2))

# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Features

* **pipeline:** Make covariate adjustment config-driven with age and sex as defaults
  * Modified step 08 to read covariates from config.yaml as a list
  * Changed default behavior to only adjust for age and sex (removed BMI and smoking from defaults)
  * Added support for `covariates_to_adjust` list in config.yaml
  * Maintained backward compatibility with legacy boolean flags
  * Added validation to ensure only valid covariates (age, sex, bmi, smoking) are used
  * Improved logging to show which covariates are being adjusted

### Bug Fixes

* **pipeline:** Fix step 05c skip logic to prevent execution when disabled
  * Fixed issue where step 05c_provenance_test executed setup code even when disabled
  * Changed from `quit()` to `stop("STEP_SKIPPED")` to properly halt execution when sourced
  * Pipeline now correctly skips step 05c and proceeds to step 05d when test_case.enabled: false

* **pipeline:** Make Youden J detection failure-proof and remove hardcoded batch references
  * Enhanced log file detection to work in both single-batch and multi-batch modes
  * Added multi-fallback approach for log file path detection
  * Removed hardcoded batch references, now uses config-driven batch determination

## [1.0.0] (2026-01-23)

### Bug Fixes

* docker tag error fixed ([6519380](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline/commit/6519380a0d4bd34e386a848fb6f1e4af60c7c03d))
* **pipeline:** make Youden J detection failure-proof and remove hardcoded batch references ([93bf990](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline/commit/93bf99015082a2297e79fb78aa26da20ae33234f))
* **release:** Fix the GitHub Actions npm cache ([c22505f](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline/commit/c22505f50f1a3609d9577e805f8028f73d75b5f7))
* Semantic Release tag detection — enhanced fix ([41c64ff](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline/commit/41c64ff34d748280945fc1a86e90585ecdfe5e1b))
* Semantic Release tag error fixed ([9cbf14b](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline/commit/9cbf14b2002a2dfd171e4d4b735d652c61b5b1db))

### Features

* **release:** implement automated release pipeline with semantic-release ([e84360c](https://github.com/FINNGEN/FinnGen3-proteomics-pipeline/commit/e84360cc11c47a16a880e3d21019af3c705ae76d))
