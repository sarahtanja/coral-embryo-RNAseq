# Bioinformatics Pipeline Roadmap
Sarah Tanja
October 7, 2024

- [Sequence files background](#sequence-files-background)
- [Coding resources](#coding-resources)
- [Handling Large Files with Git](#handling-large-files-with-git)
- [Pipeline birds-eye view](#pipeline-birds-eye-view)
  - [1. Receive raw FASTA files](#1-receive-raw-fasta-files)
  - [2. QAQC FASTA files](#2-qaqc-fasta-files)
  - [3. Align to reference genome](#3-align-to-reference-genome)
  - [4. Assemble](#4-assemble)
  - [5. Create gene expression count
    matrix](#5-create-gene-expression-count-matrix)
  - [6. Identify differentially expressed genes
    (DEGs)](#6-identify-differentially-expressed-genes-degs)
  - [7.](#7)

# Sequence files background

We’ve extracted totalRNA and shipped it to Azenta for sequencing… so
what happens next? How did we go from working with microcentrifuge tubes
in the lab to large files on the computer?

First they run the total RNA through quality assessment and quality
control steps (QAQC) to make sure it is sufficient in quantity and
quality to sequence. Here is the QAQC report from Azenta:

The total RNA is then subjected to library prep, where the RNA is turned
into cDNA.

The cDNA is what is actually sequenced, with an Illumina sequencer, 20
million reads, Poly-A selection,

The raw FASTA files come back `demultiplexed`.

# Coding resources

This roadmap was built off the following resources and references from:

- Sam White

  - [Sam’s Notebook](https://robertslab.github.io/sams-notebook/)
    - [QAQC with FastQC, MultiQC &
      Fastp](https://robertslab.github.io/sams-notebook/posts/2024/2024-10-05-FastQC-Trimming-and-QC---A.pulchra-RNA-seq-from-Azenta-Project-30-1047560508-Using-fastp/)
    - [Genome alignment with
      HISAT2](https://robertslab.github.io/sams-notebook/posts/2024/2024-10-08-RNA-seq-Alignment---A.pulchra-RNA-seq-Alignments-Using-HISAT2-and-StringTie-for-Azenta-Project-30-1047560508/index.html)

- Babraham Bioinformatics Institute:

  - [Training
    Courses](https://www.bioinformatics.babraham.ac.uk/training.html)

- Erin Chille:

  - [Mcapitata_OA_Developmental_Gene_Expression_Timeseries](https://github.com/echille/Mcapitata_OA_Developmental_Gene_Expression_Timeseries)
    github repository from the \[@chille2022\]paper

- Steven Roberts:

  - [Bioinformatics FISH
    546](https://sr320.github.io/course-fish546-2023/schedule.html)
    course at UW
  - [Lab Handbook](https://robertslab.github.io/resources/)

- Savanah Liedholt:

  - [Differentially Expressed Genes from a Reference
    Genome](https://notion.so/Differentially-Expressed-Genes-From-Reference-Genome-35267b46f51d4383b4b70bb3796d28c9)
  - [SmallRNA processing and pipeline for viral community
    Analysis](https://www.notion.so/SmallRNA-processing-and-pipeline-for-viral-community-Analysis-1f7f48c35d9e482597685354222c1852)

- Ariana Huffmyer:

  - [EarlyLifeHistory_Energetics
    TagSeq](https://github.com/AHuffmyer/EarlyLifeHistory_Energetics/tree/master)
    github repository

- Sarah Tanja:

  - [FISH 546 Compendium](https://rpubs.com/sarah_tanja/1048844)

# Handling Large Files with Git

Execute this code in terminal prior to committing to add all files under
50M

> [!TIP]
>
> `find . -type f -size -50M -exec git add {} +`

- `find .`: Starts the search from the current directory.

- `-type f`: Finds only files (not directories).

- `-size -50M`: Limits the size to files smaller than 50 megabytes.

- `-exec git add {} +`: Executes the `git add` command on all found
  files, adding them to the Git staging area.

# Pipeline birds-eye view

## 1. Receive raw FASTA files

- files are already `demultiplexed`

- files have a `.fasta.gz` zipped format

- files must be checked to make sure there were no errors in the
  transfer process (this is done with `md5sum` )

## 2. QAQC FASTA files

> [!TIP]
>
> Some great examples of previous QAQC scripts generated for *Montipora
> capitata* RNA-seq data by [E. Chille (QAQC
> Script)](https://github.com/echille/Mcapitata_OA_Developmental_Gene_Expression_Timeseries/blob/main/2-QC-Align-Assemble/mcap_rnaseq_analysis.md),
> [Sam White (Notebook
> post)](https://robertslab.github.io/sams-notebook/posts/2024/2024-10-05-FastQC-Trimming-and-QC---A.pulchra-RNA-seq-from-Azenta-Project-30-1047560508-Using-fastp/)
> and [A. Huffmyer (QAQC
> Script)](https://github.com/AHuffmyer/EarlyLifeHistory_Energetics/blob/master/Mcap2020/Scripts/TagSeq/Genome_V3/TagSeq_BioInf_genomeV3.md)

- Quality check raw sequences with
  [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) ,
  synthesize a report with [MultiQC](https://multiqc.info/)

- Clean up sequences with [Fastp](https://github.com/OpenGene/fastp)

  - Trim sequence lengths

  - Filter out bad quality reads

  - Remove adapters & polyA tails

- Check cleaned sequences with FastQC, synthesize a report with MultiQC

- Repeat cleaning steps if needed

## 3. Align to reference genome

> [!NOTE]
>
> Great example of previous HISAT2 RNA-seq alignment in Sam White’s
> notebook
> [here](https://robertslab.github.io/sams-notebook/posts/2024/2024-10-08-RNA-seq-Alignment---A.pulchra-RNA-seq-Alignments-Using-HISAT2-and-StringTie-for-Azenta-Project-30-1047560508/index.html)

- Obtain reference genome assembly & GFF annotation file
  - Genome [Version 3](http://cyanophora.rutgers.edu/montipora/) from
    Rutgers \[@stephens2022\]
  - [GFF](http://cyanophora.rutgers.edu/montipora/Montipora_capitata_HIv3.genes.gff3.gz)
    from Rutgers (or [GFF
    fixed](https://github.com/AHuffmyer/EarlyLifeHistory_Energetics/raw/master/Mcap2020/Data/TagSeq/Montipora_capitata_HIv3.genes_fixed.gff3.gz)
    from AHuffmyer?)
  - [genomes, indexes, & feature
    tracks](https://robertslab.github.io/resources/Genomic-Resources/#montipora-capitata)
    from Roberts Lab Handbook
- Align sequences to genome
  - using [HISAT-2](https://daehwankimlab.github.io/hisat2/hisat-3n/)
    \[@zhang_rapid_2021\]

## 4. Assemble

- `StringTie 2`

## 5. Create gene expression count matrix

- `Trinity`

## 6. Identify differentially expressed genes (DEGs)

- `DESeq-2`

## 7.
