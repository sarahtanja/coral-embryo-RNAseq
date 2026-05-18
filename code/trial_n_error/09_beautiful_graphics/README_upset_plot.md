# Upset Plot: Overlap of Significant DEGs

This directory contains scripts to generate an upset plot visualizing the overlap of significant differentially expressed genes (DEGs) identified across different DESeq2 analyses.

## Files

- `10_upset_plot.qmd` - Quarto document with full analysis and documentation
- `generate_upset_plot.R` - Standalone R script for quick regeneration of the plot

## Input Data

The analysis uses three sets of significant genes from DESeq2 analyses:

1. **LRT (Likelihood Ratio Test)** results: `output/06_deg/interactions/LRT_sig_genes.csv`
   - Tests for genes with differential expression across the full model
   - 44 significant genes

2. **Interaction effects**: `output/06_deg/interactions/sig_interaction_genes.csv`
   - Genes showing significant interaction between developmental stage and leachate concentration
   - 44 significant genes

3. **Main leachate effects**: `output/09_deg_multifactorial/sig_leachate.csv`
   - Genes showing significant main effect of leachate concentration
   - 2 significant genes

## Output

- **Upset plot**: `output/09_deg_multifactorial/upset_plot_deg_overlap.png`
  - 3000 x 1800 pixel PNG image at 300 DPI
  - Shows the intersection sizes and set sizes for all three categories

## Results Summary

- Total unique genes across all analyses: 44
- LRT & Interaction overlap: 42 genes
- All three categories (LRT & Interaction & Leachate): 2 genes

This shows that most genes with significant interactions are also significant in the LRT test, and a small subset also shows main leachate effects.

## Usage

### Using the Quarto document

```bash
cd code/09_beautiful_graphics
quarto render 10_upset_plot.qmd
```

### Using the standalone R script

```bash
cd code/09_beautiful_graphics
Rscript generate_upset_plot.R
```

## Dependencies

- R (>= 4.3.3)
- tidyverse
- UpSetR

Install dependencies in R:
```r
install.packages(c("tidyverse", "UpSetR"))
```

Or via apt on Ubuntu/Debian:
```bash
sudo apt-get install r-cran-tidyverse r-cran-upsetr
```
