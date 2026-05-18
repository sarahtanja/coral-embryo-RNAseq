# Gene Cluster GO Enrichment Analysis

## Overview

The file `08_cluster_GO.qmd` performs functional enrichment analysis on the 3 gene expression clusters identified in the clustering analysis (`07_cluster.qmd`).

## Purpose

This analysis characterizes the biological processes associated with each of the 3 distinct gene expression patterns observed in response to leachate exposure. By performing GO enrichment separately for each cluster, we can understand what biological functions are associated with each response pattern.

## Input Files

1. **sig_all.csv** (151 significant genes)
   - Location: `../../output/06_deg/results/sig_all.csv`
   - Contains all genes significantly affected by leachate

2. **gene_clusters.csv** (149 clustered genes)
   - Location: `../../output/07_cluster/gene_clusters.csv`
   - Maps each gene to one of 3 expression clusters
   - Cluster 1: 97 genes
   - Cluster 2: 28 genes
   - Cluster 3: 24 genes

3. **all_genes.csv**
   - Location: `../../output/06_deg/interactions/all_genes.csv`
   - Background gene set for enrichment analysis

4. **EggNog functional annotation**
   - Location: `../../input/annotations/Montipora_capitata_HIv3.genes.EggNog_results.txt`
   - Provides GO term annotations for genes

5. **Genome files** (optional for initial outline)
   - `Montipora_capitata_HIv3.assembly.fasta` - for sequence extraction
   - `Montipora_capitata_HIv3.genes_fixed.gff3` - for gene lengths (bias correction)

## Analysis Workflow

### 1. Data Loading and Annotation
- Load significant genes and cluster assignments
- Match genes to functional annotations from EggNog
- Filter for genes with GO term information

### 2. Cluster-Specific Subsets
- Separate annotated genes into 3 cluster-specific groups
- Create gene vectors for enrichment testing

### 3. GO Enrichment Analysis
- Perform GO enrichment for each cluster using `goseq`
- Apply Benjamini-Hochberg FDR correction
- Filter for significantly enriched terms (BH-adjusted p < 0.05)
- Focus on Biological Process (BP) terms

### 4. Redundancy Reduction
- Calculate semantic similarity between GO terms
- Reduce redundant terms using `rrvgo`
- Identify representative terms for each cluster

## Output Files

All outputs are saved to `../../output/08_go/`:

### Main Results
- `sig_all_clusters.csv` - Annotated genes with cluster assignments
- `GO.BP.cluster1.csv` - Enriched BP terms for cluster 1
- `GO.BP.cluster2.csv` - Enriched BP terms for cluster 2
- `GO.BP.cluster3.csv` - Enriched BP terms for cluster 3

### GO Term Mappings
- `GO.terms.cluster1.csv` - Gene-to-GO mappings for cluster 1
- `GO.terms.cluster2.csv` - Gene-to-GO mappings for cluster 2
- `GO.terms.cluster3.csv` - Gene-to-GO mappings for cluster 3

### Reduced Terms (redundancy removed)
- `reducedTerms_cluster1.csv`
- `reducedTerms_cluster2.csv`
- `reducedTerms_cluster3.csv`

### Similarity Matrices
- `simMatrix_cluster1.csv`
- `simMatrix_cluster2.csv`
- `simMatrix_cluster3.csv`

## Key Differences from 08_GO.qmd

1. **Focus**: Analyzes 3 expression clusters vs. interaction patterns
2. **Input**: Uses `gene_clusters.csv` from clustering analysis
3. **Gene set**: Works with all 151 significant genes (not just interaction genes)
4. **Flexibility**: Genome files are optional - can run without length bias correction

## Running the Analysis

### Prerequisites
```r
# Required R packages
library(GenomicRanges)
library(Biostrings)
library(goseq)
library(rrvgo)
library(org.Ce.eg.db)
library(tidyverse)
```

### Execution
```bash
# In R or RStudio
quarto render code/08_go/08_cluster_GO.qmd
```

Or open in RStudio and knit/render the document.

## Notes

- The analysis uses C. elegans database (`org.Ce.eg.db`) for GO term semantics
- Length bias correction is recommended but not required for initial analysis
- If genome files are not available, the analysis will proceed without length bias correction
- Enrichment is performed using the Wallenius approximation method in `goseq`

## Next Steps

After running this analysis:
1. Compare enriched processes across the 3 clusters
2. Visualize results using plots (e.g., in `09_beautiful_graphics/`)
3. Integrate findings with developmental stage effects
4. Interpret biological significance of each expression pattern

## Related Files

- **07_cluster.qmd** - Upstream clustering analysis that generates `gene_clusters.csv`
- **08_GO.qmd** - Parallel GO analysis for interaction patterns
- **09_beautiful_graphics/** - Downstream visualization scripts
