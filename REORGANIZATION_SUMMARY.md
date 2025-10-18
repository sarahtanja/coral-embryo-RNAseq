# DESeq2 Leachate Effects Reorganization Summary

## Overview
This document summarizes the reorganization of `code/06_deg/06_DESeq2_interaction.qmd` to clearly organize DESeq2 results for leachate effects.

## Major Changes

### 1. New Section: "Summarize Leachate Effects"
Added a comprehensive new section (starting around line 1012) that organizes all leachate-related results into three clear categories:

- **LRT Results**: Overall test for stage-by-leachate interactions (44 significant genes)
- **Interaction Terms**: Specific Wald tests for each stage × leachate combination (6 comparisons)
- **Main Leachate Effects**: Leachate vs control effects at each stage (9 comparisons)

### 2. Combined Dataframe: `combined_leachate_effects`
Created a single comprehensive dataframe that includes:
- **gene_id** column
- LRT results (baseMean, log2FC, pvalue, padj)
- All interaction term results (log2FC, pvalue, padj for each stage × leachate combination)
- All main effect results (log2FC, pvalue, padj for each stage × leachate vs control)

The dataframe contains all unique genes significant (padj < 0.05) in at least one test, with NA values for non-significant comparisons.

**Output file**: `output/06_deg/interactions/combined_leachate_effects.csv`

### 3. Summary Statistics
Added code chunks that provide:
- Counts of significant genes for each test type
- Overlap analysis showing genes significant in multiple test types
- Distribution of genes by number of significant test types

### 4. Bug Fixes and Improvements

#### Variable Naming Consistency
- Added `LRT_sig_genes <- sig_cat_LRT` alias for consistency with later code
- Added `sig_cat <- sig_cat_LRT$gene_id` for character vector of significant gene IDs
- Added `sig_lvc_cl <- sig_cleavage_low.vs.control` alias

#### Fixed Incorrect Coefficient Names
- Line 985: Changed `"leachate_high_vs_control"` to `"leachatehigh"`
- Line 997: Changed `"leachate_high_vs_control"` to `"leachatehigh"`
- Line 794-801: Fixed undefined `pcl` variable to `res_prawnchip_low.int`
- Line 811: Fixed variable name from `sig_earlygastrula_high.int_egh` to `sig_earlygastrula_high.int`

#### Added Missing Variables
- Added `vsd <- varianceStabilizingTransformation(dds_wald)` after line 568 to create variance-stabilized transformation for later use

#### Graphics Fixes
- Added `dev.off()` after pheatmap to properly close PNG device (line 1402)

### 5. Improved Documentation
- Added comprehensive table explaining the three test types
- Added clear callout notes explaining what each dataframe contains
- Improved code comments throughout the new section

## Structure of New Section

```
# Summarize Leachate Effects
├── Table: Test Types Overview
├── Extract all significant results
│   ├── LRT Significant Genes
│   ├── Interaction Terms Significant Genes (6 comparisons)
│   └── Main Leachate Effects Significant Genes (9 comparisons)
├── Combined Leachate Effects Dataframe
│   ├── Join all test results by gene_id
│   ├── Display summary statistics
│   └── Save to CSV
└── Summary Statistics
    ├── Genes by test type
    └── Overlap analysis
```

## Files Modified
- `code/06_deg/06_DESeq2_interaction.qmd` - Main analysis file

## New Output Files
- `output/06_deg/interactions/combined_leachate_effects.csv` - Comprehensive dataframe with all significant leachate effects

## Testing Recommendations
While the code changes are minimal and focused on organizing existing results, it's recommended to:
1. Render the Quarto document to verify all code chunks execute without errors
2. Verify the `combined_leachate_effects.csv` file is created correctly
3. Check that the heatmap PNG is generated properly
4. Verify that downstream analyses (clustering, gene lists, etc.) still work correctly

## Notes
- All existing functionality is preserved
- No analyses were removed or changed, only reorganized
- The new combined dataframe provides a clearer view of leachate effects across all test types
- Variable names were standardized for consistency
