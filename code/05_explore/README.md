# Step 5: Data Exploration Figures

This folder contains the initial exploratory data analysis of the gene count matrix for the *Montipora capitata* embryo RNA-seq experiment examining the effects of PVC leachate exposure.

## Overview

The exploratory analysis includes:
- Quality control checks on the count matrix
- Principal Component Analysis (PCA) to visualize sample clustering
- Heatmap analysis to examine gene expression patterns
- Statistical tests for batch effects (betadisper and PERMANOVA)

See [`05_data_exploration.qmd`](05_data_exploration.qmd) for the full analysis code and [`05_data_exploration.md`](05_data_exploration.md) or [`05_data_exploration.html`](05_data_exploration.html) for the rendered output.

---

## Main Output Figures

### Principal Component Analysis (PCA)

![PCA Plot](../../output/05_explore/pca.png)

*Principal component analysis showing sample clustering by treatment and developmental stage.*

### Z-Score Heatmap

![Z-Score Heatmap](../../output/05_explore/zscore_heatmap.png)

*Heatmap showing hierarchical clustering of samples based on gene expression patterns (z-score transformed).*

### PCA & Heatmap Combined

![PCA & Heatmap Full Gene](../../output/05_explore/fig_pca_heatmap_fullgene.png)

*Combined visualization of PCA and heatmap for the full gene set.*

### Z-Score Heatmap (Alternative)

![Z-Score Heatmap](images/zscore_heatmap.png)

*Alternative z-score heatmap visualization.*

---

## Generated Figures from Analysis

The following figures are automatically generated from the Quarto analysis document and show various quality control, filtering, and exploratory analysis steps.

### Distribution of Unfiltered Reads

![Distribution of Unfiltered Reads](05_data_exploration_files/figure-commonmark/distr_unfiltered_reads-1.png)

### Quality Control and Filtering Steps

![QC Step 1](05_data_exploration_files/figure-commonmark/unnamed-chunk-8-1.png)

![QC Step 2](05_data_exploration_files/figure-commonmark/unnamed-chunk-10-1.png)

![QC Step 3](05_data_exploration_files/figure-commonmark/unnamed-chunk-12-1.png)

![QC Step 4](05_data_exploration_files/figure-commonmark/unnamed-chunk-13-1.png)

![QC Step 5](05_data_exploration_files/figure-commonmark/unnamed-chunk-14-1.png)

![QC Step 6](05_data_exploration_files/figure-commonmark/unnamed-chunk-15-1.png)

![QC Step 7](05_data_exploration_files/figure-commonmark/unnamed-chunk-17-1.png)

![QC Step 8](05_data_exploration_files/figure-commonmark/unnamed-chunk-18-1.png)

![QC Step 9](05_data_exploration_files/figure-commonmark/unnamed-chunk-19-1.png)

![QC Step 10](05_data_exploration_files/figure-commonmark/unnamed-chunk-20-1.png)

![QC Step 11](05_data_exploration_files/figure-commonmark/unnamed-chunk-22-1.png)

![QC Step 12](05_data_exploration_files/figure-commonmark/unnamed-chunk-23-1.png)

![QC Step 13](05_data_exploration_files/figure-commonmark/unnamed-chunk-24-1.png)

![QC Step 14](05_data_exploration_files/figure-commonmark/unnamed-chunk-25-1.png)

![QC Step 15](05_data_exploration_files/figure-commonmark/unnamed-chunk-27-1.png)

![QC Step 16](05_data_exploration_files/figure-commonmark/unnamed-chunk-28-1.png)

![QC Step 17](05_data_exploration_files/figure-commonmark/unnamed-chunk-31-1.png)

![QC Step 18](05_data_exploration_files/figure-commonmark/unnamed-chunk-32-1.png)

![QC Step 19](05_data_exploration_files/figure-commonmark/unnamed-chunk-34-1.png)

![QC Step 20](05_data_exploration_files/figure-commonmark/unnamed-chunk-36-1.png)

![QC Step 21](05_data_exploration_files/figure-commonmark/unnamed-chunk-41-1.png)

![QC Step 22](05_data_exploration_files/figure-commonmark/unnamed-chunk-43-1.png)

![QC Step 23](05_data_exploration_files/figure-commonmark/unnamed-chunk-45-1.png)

![QC Step 24](05_data_exploration_files/figure-commonmark/unnamed-chunk-47-1.png)

![QC Step 25](05_data_exploration_files/figure-commonmark/unnamed-chunk-51-1.png)

![QC Step 26](05_data_exploration_files/figure-commonmark/unnamed-chunk-52-1.png)

---

## Inputs

- `output/04_count/gene_count_matrix.csv` - Unfiltered gene count matrix
- `metadata/metadata.csv` - Sample metadata

## Outputs

- `output/05_explore/gcm_filtor.csv` - Filtered gene count matrix (gcm) with outlier samples removed
- `metadata/metadata_or.csv` - Filtered metadata with outlier samples removed
- `output/05_explore/pca.png` - Exploratory PCA plot
- `output/05_explore/zscore_heatmap.png` - Exploratory heatmap
- `output/05_explore/fig_pca_heatmap_fullgene.png` - Combined PCA and heatmap visualization
- `output/05_explore/permanova_result_fullgeneset.csv` - PERMANOVA test results
- `output/05_explore/vst_matrix.rds` - Variance stabilized transformation matrix
- `output/05_explore/pheat_annotation_df.rds` - Heatmap annotation data frame

---

## Analysis Summary

This exploratory analysis:
1. Filters the gene count matrix to remove low gene counts
2. Creates a DESeq2 Dataset Object for visualization
3. Identifies and removes outlier samples through PCA
4. Examines potential batch effects (night, leachate, stage)
5. Performs PERMANOVA to test for significant differences
6. Generates heatmaps to visualize broad expression patterns

For detailed methodology and results, see the [analysis document](05_data_exploration.md).
