# FinnGen 3 Batch 2 Olink Proteomics Data Release Note

- **Release Date**: January 2026 (v.2.1)
- **Platform**: Olink Explore HT (5K)
- **Batch**: FG3 Batch 2
- **Total Measurements (raw import)**: 14,144,000 (2,600 samples × 5,440 proteins)
- **Total Measurements (after QC)**: 13,282,832 (2,452 samples × 5,416 proteins)
- **Author**: Reza Jabal, PhD (rjabal@broadinstitute.org)
- **Reviewer**: Mitja Kurki, PhD (mkurki@broadinstitute.org)

# -----------------------------------------------------------

## File Descriptions

### Comprehensive Outliers List

**Files**:
- `comprehensive_outliers_list_fg3_batch_02.tsv` (tab-separated format)
- `comprehensive_outliers_list_fg3_batch_02.parquet` (Parquet format)

**Description**: Contains all **75 samples** flagged as outliers by any QC method, with detailed QC flags and metrics.

### Clean Proteomics Dataset

**Files**: `npx_matrix_all_qc_passed_fg3_batch_02.rds`, `.parquet`, and `.tsv`

**Dimensions**: **2,452 samples** (rows) × 5,416 proteins (columns)

## Sample QC Summary

### Analysis-Ready Dataset
- **Analysis-ready samples**: 2,527 (2,477 FinnGen + 50 Bridge, after excluding 47 non-FinnGen samples)

## Technical Outlier Detection

### PCA Outlier Detection
- **Input samples**: 2,527 (analysis-ready samples)
- **Constant/zero-variance proteins removed**: 0 (for PCA calculation only)
- **Sequential filtering results**:
  - PC1/PC2 filtering (±5 SD): 0 samples flagged
  - PC3/PC4 filtering (±5 SD): 17 samples flagged
  - Median filtering (±5 SD): 0 additional samples flagged
  - IQR filtering (±5 SD): 5 samples flagged
- **Total PCA outliers flagged**: **22 samples** (0.87% of 2,527 analysis-ready samples)

### Technical Outlier Detection
- **Input samples**: 2,527 (analysis-ready samples)
- **Total technical outliers flagged**: **27 samples** (1.07% of 2,527 analysis-ready samples)
  - SD outliers: 22 samples
  - Missing rate (>5%): 5 samples

### Z-score Outlier Detection
- **Input samples**: 2,527 (analysis-ready samples)
- **Z-score outliers flagged**: **7 samples** (0.28% of 2,527 analysis-ready samples)

## Provenance Steps

### Sex Outlier Detection
- **Input samples**: 2,505 (2,527 − 22 PCA outliers removed)
- **Strict mismatches flagged**: **17 samples** (0.68% of 2,505 PCA-cleaned samples)
- **Sex outliers flagged**: **5 samples** (0.20% of 2,505 PCA-cleaned samples)
- **Total sex-related flags**: **22 samples** (0.88% of 2,505 PCA-cleaned samples)

### pQTL-based Outlier Detection
- **Input samples**: 2,505 (2,527 − 22 PCA outliers removed)
- **Outliers flagged (MeanAbsZ-based)**: **14 samples** (0.56% of 2,505 PCA-cleaned samples)

### Final QC Integration
- **Unique samples flagged**: **75 samples** (2.97% of 2,527 analysis-ready samples)
- **Samples flagged by multiple methods**: **15 samples** (20.0% of all flagged samples)
  - 2 methods: 8 samples
  - 3 methods: 7 samples

**Breakdown by QC method** (all percentages relative to 2,527 analysis-ready samples):
- Initial QC: 5 samples (0.20%)
- PCA: 22 samples (0.87%)
- Sex Mismatch (Strict): 17 samples (0.67%)
- Sex Outlier (Threshold): 5 samples (0.20%)
- Technical: 27 samples (1.07%)
- Z-score: 7 samples (0.28%)
- pQTL: 14 samples (0.55%)

**Final clean dataset**:
- **Samples passing all QC**: **2,452 samples** (97.03% retention rate from 2,527 analysis-ready samples)
- **Biological proteins**: 5,416 (24 control probes excluded from released data)
- **Clean NPX matrix**: 2,452 samples × 5,416 proteins

## Overall Sample Flow Summary

```
Raw parquet:              2,600 samples (100%)
  ↓ Filter controls       -20 (removed)
Biological samples:       2,580 samples (100%)
  ↓ Initial QC            -6 (5 FINNGEN + 1 non-FinnGen, removed)
After QC:                 2,574 samples (99.8%)
  ↓ Exclude non-FinnGen   -47 (excluded from analysis-ready)
Analysis-ready:           2,527 samples (98.0%)
  ↓ TECHNICAL OUTLIER DETECTION (Parallel flagging on base matrix)
     PCA:                  22 flagged (0.87%)
     Technical:            27 flagged (1.07%)
     Z-score:              7 flagged (0.28%)
  ↓ PCA-cleaned matrix:    2,505 samples (2,527 - 22 PCA outliers)
  ↓ PROVENANCE STEPS (Sequential on PCA-cleaned matrix)
     Sex:                  22 flagged (0.88%: 17 strict + 5 threshold)
     pQTL:                 14 flagged (0.56%)
  ↓ Final QC Integration: Combine flags (union logic)
     Unique samples flagged: 75 (2.97% of 2,527)
     Overlaps: 15 samples flagged by multiple methods
  ↓ Final QC Integration: Remove all flagged samples
Final (pre-normalisation): 2,452 samples (97.03% of 2,527 analysis-ready)
```

**Retention rate**: 97.03% (2,452/2,527 analysis-ready samples) or 95.04% (2,452/2,580 biological samples)

## Quality Control Thresholds Summary

| Method | Threshold | Rationale |
|--------|-----------|-----------|
| **PCA (PC1/PC2)** | mean ± 5 × SD (after Olink scaling) | Matches Batch 1, ~99.9999% specificity |
| **PCA (PC3/PC4)** | mean ± 5 × SD | Sequential filtering, union of flags |
| **PCA (Sample median)** | mean ± 5 × SD | Detects extreme central tendency |
| **PCA (Sample IQR)** | mean ± 5 × SD | Detects unusual variability patterns |
| **Sex mismatch** | predicted_sex ≠ genetic_sex (0.5 threshold) | Binary label error, severe classification errors |
| **Sex outlier** | 0.5 threshold | Borderline predictions, mild warning |
| **Technical (Plate/Batch/Processing)** | 5 × MAD ≈ 4 × SD | Robust, harmonised with z-score method |
| **Technical (Sample mean NPX)** | 5 × MAD ≈ 4 × SD | Two-sided, robust to outliers |
| **Technical (Sample SD NPX)** | median + 4 × MAD | One-sided upper, high variance detection |
| **Technical (Missing rate)** | 5% (fixed) | More stringent than Initial QC, catches borderline cases |
| **Z-score (Per-protein)** | \|Z\| > 4 | ~99.994% specificity, harmonised with Technical method |
| **Z-score (Sample threshold)** | >10% proteins | Requires systematic issues, not isolated extremes |
| **pQTL (MeanAbsZ)** | mean + 4 × SD | Population-based, SD-based (not MAD), exclusively used for outlier assignment |

## Notes

**Document version v.2.1 update (January 21, 2026)**: This version reflects the refactored pipeline implementation after initial QC integration fix. Key differences from v.02:
- Analysis-ready samples: 2,527 (vs 2,522 in v.02) - 5 additional samples
- Initial QC: 5 samples tracked (vs 0 in v.02) - now correctly integrated
- PCA outliers: 22 (vs 24 in v.02) - 2 fewer outliers
- Technical outliers: 27 (vs 27 in v.02) - matches expected
- Sex mismatch outliers: 17 (vs 18 in v.02) - 1 fewer outlier
- Sex threshold outliers: 5 (vs 13 in v.02) - 8 fewer outliers (investigation ongoing)
- Total outliers: 75 (vs 86 in v.02) - 11 fewer outliers overall
- Final QC passed samples: 2,452 (vs 2,441 in v.02) - 11 additional samples
- Multiple method overlaps: 15 (vs 15 in v.02) - matches expected

These differences are primarily due to:
- Initial QC integration fix: 5 outliers now correctly tracked
- Different input sample counts (2,527 vs 2,522) affecting threshold calculations
- Parallel flagging implementation enabling consistent detection
- Different constant protein removal (0 vs 8) affecting PCA calculations
- Sex outlier detection discrepancy: 5 found vs 13 expected (investigation ongoing - may be due to Youden threshold difference: 0.63 vs 0.71)

All thresholds match the original implementation (5×SD for PCA, 5×MAD for technical, 4×SD for Z-score and pQTL). Pipeline version: v1.5.0

---

**Pipeline Version**: v1.5.0 (January 2026 - Refactored)
**Data Release**: FG3 Batch 2
**Final Sample Count**: 2,452 samples × 5,416 proteins (24 control probes excluded from released data)
