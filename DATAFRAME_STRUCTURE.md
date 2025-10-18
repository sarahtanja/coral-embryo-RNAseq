# Combined Leachate Effects Dataframe Structure

The `combined_leachate_effects.csv` dataframe has the following structure:

## Columns

### Gene Identifier
- `gene_id` - Unique gene identifier

### LRT Results (4 columns)
- `LRT_baseMean` - Mean normalized counts across all samples
- `LRT_log2FC` - Log2 fold change from LRT
- `LRT_pvalue` - Raw p-value from LRT
- `LRT_padj` - BH-adjusted p-value from LRT

### Interaction Terms (18 columns = 6 comparisons × 3 metrics)
Each comparison has log2FC, pvalue, and padj:
- `int_prawnchip_low_*` - Prawnchip × Low interaction
- `int_prawnchip_mid_*` - Prawnchip × Mid interaction
- `int_prawnchip_high_*` - Prawnchip × High interaction
- `int_earlygastrula_low_*` - Early Gastrula × Low interaction
- `int_earlygastrula_mid_*` - Early Gastrula × Mid interaction
- `int_earlygastrula_high_*` - Early Gastrula × High interaction

### Main Effects (27 columns = 9 comparisons × 3 metrics)
Each comparison has log2FC, pvalue, and padj:

**Cleavage Stage:**
- `main_cleavage_low_*` - Cleavage Low vs Control
- `main_cleavage_mid_*` - Cleavage Mid vs Control
- `main_cleavage_high_*` - Cleavage High vs Control

**Prawnchip Stage:**
- `main_prawnchip_low_*` - Prawnchip Low vs Control
- `main_prawnchip_mid_*` - Prawnchip Mid vs Control
- `main_prawnchip_high_*` - Prawnchip High vs Control

**Early Gastrula Stage:**
- `main_earlygastrula_low_*` - Early Gastrula Low vs Control
- `main_earlygastrula_mid_*` - Early Gastrula Mid vs Control
- `main_earlygastrula_high_*` - Early Gastrula High vs Control

## Total Columns
**50 columns total**: 1 (gene_id) + 4 (LRT) + 18 (interactions) + 27 (main effects)

## Rows
Each row represents a unique gene that is significant (BH-adjusted p-value < 0.05) in at least one test:
- LRT test, OR
- Any interaction term, OR
- Any main effect comparison

## NA Values
- NA values indicate the gene was **not significant** (padj ≥ 0.05) for that particular test
- This allows quick filtering to genes significant for specific comparisons

## Example Usage

### Find genes with significant LRT
```r
combined_leachate_effects %>%
  filter(!is.na(LRT_padj))
```

### Find genes with any significant interaction
```r
combined_leachate_effects %>%
  filter(if_any(matches("^int_.*_padj$"), ~ !is.na(.)))
```

### Find genes with cleavage stage effects
```r
combined_leachate_effects %>%
  filter(if_any(matches("^main_cleavage_.*_padj$"), ~ !is.na(.)))
```

### Find genes significant in multiple test types
```r
combined_leachate_effects %>%
  mutate(
    has_LRT = !is.na(LRT_padj),
    has_interaction = rowSums(!is.na(select(., matches("^int_.*_padj$")))) > 0,
    has_main_effect = rowSums(!is.na(select(., matches("^main_.*_padj$")))) > 0,
    n_test_types = has_LRT + has_interaction + has_main_effect
  ) %>%
  filter(n_test_types > 1)
```
