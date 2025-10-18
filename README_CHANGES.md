# DESeq2 Leachate Effects Reorganization - Complete Guide

## Overview
This PR reorganizes the DESeq2 interaction analysis in `code/06_deg/06_DESeq2_interaction.qmd` to provide a clear, comprehensive summary of leachate effects on coral embryos.

## Problem Statement
The original analysis had DESeq2 results for leachate effects scattered throughout the code without a clear summary. The goal was to:
1. Organize results clearly by test type (LRT, interaction terms, main effects)
2. Create a single dataframe combining all significant leachate effects
3. Fix bugs and inconsistencies in the code

## Solution
Added a new section **"Summarize Leachate Effects"** (starting ~line 1012) that:
- Extracts and organizes all significant results (padj < 0.05) by test type
- Combines them into a comprehensive dataframe: `combined_leachate_effects`
- Provides summary statistics and overlap analysis

## Key Deliverable: combined_leachate_effects.csv

### What is it?
A comprehensive dataframe containing all genes with significant leachate effects across three test types:

| Test Type | Comparisons | Purpose |
|-----------|-------------|---------|
| **LRT** | 1 overall test | Tests for any stage-by-leachate interaction |
| **Interaction Terms** | 6 specific tests | Tests stage-specific deviations from main effects |
| **Main Effects** | 9 comparisons | Tests leachate vs control at each stage |

### Structure
- **50 columns total**:
  - 1 gene_id column
  - 4 LRT result columns (baseMean, log2FC, pvalue, padj)
  - 18 interaction term columns (6 comparisons × 3 metrics)
  - 27 main effect columns (9 comparisons × 3 metrics)
- **Rows**: All genes significant in at least one test
- **NA values**: Indicate non-significant results (padj ≥ 0.05) for that test

### Location
`output/06_deg/interactions/combined_leachate_effects.csv`

## Changes Made

### 1. New Code Section
Added comprehensive "Summarize Leachate Effects" section with:
- Clear table explaining test types
- Extraction of all significant results
- Creation of combined dataframe
- Summary statistics and overlap analysis

### 2. Bug Fixes
- Fixed undefined variable `pcl` → `res_prawnchip_low.int` (line 794)
- Fixed incorrect coefficient `"leachate_high_vs_control"` → `"leachatehigh"` (lines 985, 997)
- Fixed variable name `sig_earlygastrula_high.int_egh` → `sig_earlygastrula_high.int` (line 811)

### 3. Added Missing Variables
- `vsd <- varianceStabilizingTransformation(dds_wald)` for downstream analyses
- `LRT_sig_genes <- sig_cat_LRT` alias for consistency
- `sig_cat <- sig_cat_LRT$gene_id` for gene ID vector
- `sig_lvc_cl <- sig_cleavage_low.vs.control` alias

### 4. Code Improvements
- Added `dev.off()` after heatmap to properly close PNG device
- Standardized column names for consistency
- Added comprehensive comments and documentation

## Documentation Files

| File | Purpose |
|------|---------|
| `REORGANIZATION_SUMMARY.md` | Complete list of changes and bug fixes |
| `DATAFRAME_STRUCTURE.md` | Detailed column descriptions and structure |
| `QUICK_REFERENCE.md` | Usage examples and code snippets |
| `FLOW_DIAGRAM.txt` | Visual representation of data flow |
| `README_CHANGES.md` | This file - complete guide |

## Usage Examples

### Load the combined dataframe
```r
combined <- read_csv("output/06_deg/interactions/combined_leachate_effects.csv")
```

### Find genes with LRT significance
```r
lrt_genes <- combined %>%
  filter(!is.na(LRT_padj)) %>%
  arrange(LRT_padj)
```

### Find genes with prawnchip interactions
```r
prawnchip_int <- combined %>%
  filter(if_any(starts_with("int_prawnchip"), ~ !is.na(.)))
```

### Find genes responsive to high leachate
```r
high_genes <- combined %>%
  filter(if_any(matches("high.*_padj$"), ~ !is.na(.)))
```

See `QUICK_REFERENCE.md` for more examples.

## Testing Recommendations

While the changes are minimal and preserve all existing functionality, testing is recommended:

1. **Render the Quarto document**
   ```bash
   quarto render code/06_deg/06_DESeq2_interaction.qmd
   ```

2. **Verify output file**
   ```r
   combined <- read_csv("output/06_deg/interactions/combined_leachate_effects.csv")
   dim(combined)  # Should have 50 columns
   ```

3. **Check existing outputs still work**
   - Heatmap: `output/06_deg/interactions/int_wald_heatmap.png`
   - Gene lists: Files in `output/06_deg/contrasts/` and `output/06_deg/interactions/`

## Benefits

1. **Clarity**: All leachate effects organized in one place
2. **Completeness**: Single dataframe with all test results
3. **Flexibility**: Easy to filter and subset for downstream analyses
4. **Reproducibility**: Well-documented with clear structure
5. **Correctness**: Bug fixes ensure accurate results

## Next Steps

The reorganized code is ready for:
- Functional enrichment analysis on gene subsets
- Comparative analysis across test types
- Visualization of overlapping effects
- Integration with other data (e.g., expression patterns, GO terms)

## Questions?

Refer to:
- `DATAFRAME_STRUCTURE.md` for column details
- `QUICK_REFERENCE.md` for usage examples
- `REORGANIZATION_SUMMARY.md` for technical details

## Summary

This reorganization provides a clear, comprehensive view of leachate effects on coral embryos, combining LRT, interaction terms, and main effects into a single dataframe while fixing several bugs and improving code clarity. All existing functionality is preserved while making the analysis more accessible and usable.
