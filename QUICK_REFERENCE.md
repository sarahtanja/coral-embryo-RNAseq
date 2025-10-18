# Quick Reference: Using combined_leachate_effects.csv

## What is it?
A comprehensive dataframe containing all genes with significant (BH-adjusted p-value < 0.05) leachate effects across three test types:
1. **LRT** - Likelihood Ratio Test for stage-by-leachate interactions
2. **Interaction Terms** - Wald tests for specific stage × leachate interaction effects
3. **Main Effects** - Wald contrasts for leachate vs control at each stage

## Where to find it?
`output/06_deg/interactions/combined_leachate_effects.csv`

## Key Features
- **50 columns**: 1 gene_id + 4 LRT metrics + 18 interaction metrics + 27 main effect metrics
- **NA values** indicate non-significant results (padj ≥ 0.05) for that test
- Each gene is included if significant in **at least one** test

## Common Use Cases

### 1. View all LRT-significant genes
```r
lrt_genes <- combined_leachate_effects %>%
  filter(!is.na(LRT_padj)) %>%
  arrange(LRT_padj)
```

### 2. Find genes with prawnchip interactions
```r
prawnchip_int <- combined_leachate_effects %>%
  filter(if_any(c(int_prawnchip_low_padj, 
                  int_prawnchip_mid_padj, 
                  int_prawnchip_high_padj), ~ !is.na(.)))
```

### 3. Find genes responsive to mid leachate at any stage
```r
mid_responsive <- combined_leachate_effects %>%
  filter(if_any(c(main_cleavage_mid_padj,
                  main_prawnchip_mid_padj,
                  main_earlygastrula_mid_padj), ~ !is.na(.)))
```

### 4. Identify genes significant in all three test types
```r
triple_sig <- combined_leachate_effects %>%
  mutate(
    has_LRT = !is.na(LRT_padj),
    has_int = rowSums(!is.na(select(., matches("^int_.*_padj$")))) > 0,
    has_main = rowSums(!is.na(select(., matches("^main_.*_padj$")))) > 0
  ) %>%
  filter(has_LRT & has_int & has_main)
```

### 5. Get strongest effects (most significant) for each test type
```r
# Strongest LRT effects
top_lrt <- combined_leachate_effects %>%
  filter(!is.na(LRT_padj)) %>%
  slice_min(LRT_padj, n = 10)

# Strongest interaction in early gastrula
top_eg_int <- combined_leachate_effects %>%
  filter(!is.na(int_earlygastrula_high_padj)) %>%
  slice_min(int_earlygastrula_high_padj, n = 10)
```

### 6. Create a heatmap of significant effects
```r
# Prepare data for heatmap
heatmap_data <- combined_leachate_effects %>%
  filter(!is.na(LRT_padj)) %>%
  select(gene_id, matches("_log2FC$")) %>%
  column_to_rownames("gene_id")

# Plot
pheatmap::pheatmap(heatmap_data, 
                   scale = "row",
                   cluster_rows = TRUE,
                   cluster_cols = TRUE,
                   na_value = "grey90")
```

### 7. Export subsets for downstream analysis
```r
# Export genes responsive to high leachate
high_genes <- combined_leachate_effects %>%
  filter(if_any(matches("high.*_padj$"), ~ !is.na(.))) %>%
  select(gene_id)

write_csv(high_genes, "high_leachate_responsive_genes.csv")
```

## Column Name Patterns

### LRT columns
- `LRT_baseMean`, `LRT_log2FC`, `LRT_pvalue`, `LRT_padj`

### Interaction columns
Pattern: `int_<stage>_<leachate>_<metric>`
- Stages: `prawnchip`, `earlygastrula`
- Leachate: `low`, `mid`, `high`
- Metrics: `log2FC`, `pvalue`, `padj`

### Main effect columns
Pattern: `main_<stage>_<leachate>_<metric>`
- Stages: `cleavage`, `prawnchip`, `earlygastrula`
- Leachate: `low`, `mid`, `high`
- Metrics: `log2FC`, `pvalue`, `padj`

## Tips
1. Use `matches()` or `starts_with()` from dplyr to select groups of columns
2. Use `if_any()` to filter rows with at least one significant result
3. Use `if_all()` to filter rows where all specified tests are significant
4. Remember: NA means **not significant** (padj ≥ 0.05), not missing data

## For More Information
See `DATAFRAME_STRUCTURE.md` for complete column descriptions and `REORGANIZATION_SUMMARY.md` for details on how the dataframe was created.
