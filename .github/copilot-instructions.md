# GitHub Copilot Instructions

## Repository Overview
This repository contains code and analysis for RNA-seq gene expression analysis of *Montipora capitata* (rice coral) embryos exposed to three increasing levels of polyvinyl chloride (PVC) leachate. The study examines the impact of PVC leachate on early coral development through a single-stressor ecotoxicology assay.

## Technology Stack
- **Primary Language**: R
- **Documentation**: Quarto (.qmd files)
- **Bioinformatics Tools**:
  - FastQC & Fastp: Quality control and read trimming
  - HISAT2: Read alignment to reference genome
  - StringTie: Assembly
  - DESeq2: Differential expression analysis
  - goseq: GO enrichment analysis
  - ggplot2 & rrvgo: Visualization

## Project Structure
The repository follows a sequential bioinformatics pipeline:

1. **`01_get/`**: Download and prepare reference genome and sequences
2. **`02_qaqc/`**: Quality control of raw FASTQ files
3. **`03_align/`**: Read alignment using HISAT2
4. **`04_count/`**: Generate gene expression count matrix
5. **`05_explore/`**: Data exploration and visualization
6. **`06_deg/`**: Differential expression analysis with DESeq2
7. **`07_cluster/`**: Gene clustering analysis
8. **`08_go/`**: Gene Ontology enrichment analysis
9. **`09_beautiful_graphics/`**: Publication-ready visualizations

### Data Directories
- **`rawfastq/`**: Raw sequencing data (FASTQ files)
- **`input/`**: Reference genomes and annotation files
- **`output/`**: Analysis results and intermediate files
- **`metadata/`**: Sample metadata and experimental design
- **`images/`**: Figures and plots

## Reference Genome
- **Species**: *Montipora capitata*
- **Version**: Genome Version 3 from Rutgers
- **Source**: http://cyanophora.rutgers.edu/montipora/

## Coding Standards
- All analysis scripts are written as Quarto documents (.qmd)
- Each analysis step produces HTML and markdown outputs
- Follow existing naming conventions: `NN_description.qmd` where NN is the step number
- Use relative paths from the project root
- Document code with comments explaining biological/analytical reasoning

## Key Principles
1. **Reproducibility**: All analyses should be reproducible from raw data
2. **Documentation**: Each step should have clear documentation of parameters and rationale
3. **Version Control**: Avoid committing large files (>50MB). Use Git LFS or exclude via .gitignore
4. **Sequential Processing**: Respect the pipeline order; later steps depend on earlier outputs

## Important Notes
- GO enrichment requires multiple genes; single-gene enrichment is not statistically meaningful
- For single genes in a cluster, report functional annotation directly (see RECOMMENDATIONS.md)
- Raw FASTQ files are demultiplexed and gzipped
- Use md5sum to verify file transfers

## Related Repositories
- [coral-embryo-scope](https://github.com/sarahtanja/coral-embryo-scope): Development analysis
- [coral-embryo-microbiome](https://github.com/sarahtanja/coral-embryo-microbiome): Microbiome analysis

## Helpful Resources
- Roberts Lab Handbook: https://robertslab.github.io/resources/
- Sam White's Notebook: https://robertslab.github.io/sams-notebook/
- Babraham Bioinformatics Training: https://www.bioinformatics.babraham.ac.uk/training.html
