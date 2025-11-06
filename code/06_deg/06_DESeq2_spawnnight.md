# Multifactorial Differential Gene Expression
Sarah Tanja
2024-12-02

- [<span class="toc-section-number">1</span> Background](#background)
  - [<span class="toc-section-number">1.1</span> Inputs](#inputs)
  - [<span class="toc-section-number">1.2</span> Outputs](#outputs)
- [<span class="toc-section-number">2</span> Setup](#setup)
  - [<span class="toc-section-number">2.1</span> Install
    packages](#install-packages)
  - [<span class="toc-section-number">2.2</span> Load
    packages](#load-packages)
  - [<span class="toc-section-number">2.3</span> Import filtered gene
    count matrix](#import-filtered-gene-count-matrix)
  - [<span class="toc-section-number">2.4</span> Import
    metadata](#import-metadata)
- [<span class="toc-section-number">3</span> LRT](#lrt)
- [<span class="toc-section-number">4</span> Make
  DESeqDataSet](#make-deseqdataset)
  - [<span class="toc-section-number">4.1</span> Run DESeq LRT
    Tests](#run-deseq-lrt-tests)
    - [<span class="toc-section-number">4.1.1</span> Spawn
      night](#spawn-night)
    - [<span class="toc-section-number">4.1.2</span> Interaction
      term](#interaction-term)
    - [<span class="toc-section-number">4.1.3</span> Main leachate
      effect](#main-leachate-effect)
    - [<span class="toc-section-number">4.1.4</span> Joint](#joint)
- [<span class="toc-section-number">5</span> Wald](#wald)
  - [<span class="toc-section-number">5.1</span> Results](#results)
    - [<span class="toc-section-number">5.1.1</span> Interaction
      effects](#interaction-effects)
  - [<span class="toc-section-number">5.2</span> Contrasts](#contrasts)
    - [<span class="toc-section-number">5.2.1</span>
      Cleavage](#cleavage)
    - [<span class="toc-section-number">5.2.2</span>
      Prawnchip](#prawnchip-1)
    - [<span class="toc-section-number">5.2.3</span> Early
      gastrula](#early-gastrula-1)
- [<span class="toc-section-number">6</span> Summary](#summary)
- [<span class="toc-section-number">7</span> Next steps](#next-steps)

# Background

Here we are looking at a multifactorial differential gene expression
analysis for *Montipora capitata* coral embryos at 3 developmental
stages (4 hours post fertilization = cleavage , 9 hours post
fertilization = prawnchip, and 14 hours post fertilization = early
gastrula). Embryos at each stage were exposed to either a control, low,
mid, or high concentration of PVC leachate. This work aims to understand
sublethal effects of ubiquitous non-point source PVC leachate pollution
on coral reefs.

https://chatgpt.com/share/68ed79f1-e768-8008-9caa-0df1ca68b4ef

## Inputs

- `gcm_filtor.csv` is the filtered gene count matrix (genes were kept in
  which 85% of the samples had a minimum count of 15) for 61 samples
  (two outliers removed).
- `metadata_or.csv` is the metadata for 61 samples (two outliers
  removed)

## Outputs

- `output/06_deg/resul/int_wald_heatmap.png`
- `volcano_facet_interaction.png`

# Setup

``` r
output_path <- "../../output/06_deg/spawnnight/"
metadata_path <- "../../metadata/metadata_or.csv"
```

## Install packages

``` r
# CRAN packages
#install.packages(c(
#  "tidyverse",
#  "viridis",
#  "plotly",
#  "cowplot",
#  "pheatmap",
#  "vsn",
#  "ashr",
#  "remotes"  # required for GitHub installs
#))

#install.packages("ComplexUpset")    # ggplot2-based, prettier
# GitHub package
#remotes::install_github("mdscheuerell/flexoki")

# Bioconductor
#if (!"BiocManager" %in% rownames(installed.packages())) {
#  install.packages("BiocManager")
#}
#BiocManager::install(c("DESeq2", "DEGreport"))
```

## Load packages

``` r
# Load packages
library(DESeq2)
library(viridis)
library(tidyverse)
library(DEGreport)
library(ashr)
library(flexoki)
library(plotly)
library(vsn)
library(pheatmap)
library(ComplexUpset)
library(ggplot2)
library(ggplotify)
library(patchwork)   # for wrap_plots(...)
library(cowplot)     # for plot_grid & get_legend
library(dplyr)
library(grid) 
library(purrr)
library(forcats)
library(ggforce)   # for facet_wrap_paginate
```

## Import filtered gene count matrix

In this count matrix, each row represents a gene, each column a
sequenced RNA library (i.e. a sample), and the values give the estimated
counts of fragments that were assigned to the respective gene in each
library by `HISAT2`

- `check.names = FALSE` removed the ‘X’ from in front of sample id’s in
  the column headers

``` r
gcm_filtor <- read.csv("../../output/05_explore/gcm_filtor.csv", row.names="gene_id", check.names=FALSE) 
head(gcm_filtor)
```

                                              101112C14 101112C4 101112C9 101112H14
    Montipora_capitata_HIv3___RNAseq.g4581.t1       176      458      395       195
    Montipora_capitata_HIv3___RNAseq.g4582.t1      1137      380      934       882
    Montipora_capitata_HIv3___RNAseq.g4583.t1      2043      905     1848      1847
    Montipora_capitata_HIv3___RNAseq.g4584.t1      2641      452      706      1917
    Montipora_capitata_HIv3___RNAseq.g4585.t1       521       17      161       304
    Montipora_capitata_HIv3___RNAseq.g4586.t1       121       33       81       116
                                              101112H4 101112H9 101112L14 101112L4
    Montipora_capitata_HIv3___RNAseq.g4581.t1      144      351       156      187
    Montipora_capitata_HIv3___RNAseq.g4582.t1      326     1337      1101      392
    Montipora_capitata_HIv3___RNAseq.g4583.t1      945     2087      1625     1039
    Montipora_capitata_HIv3___RNAseq.g4584.t1      266      625      2301      367
    Montipora_capitata_HIv3___RNAseq.g4585.t1       33      243       449       20
    Montipora_capitata_HIv3___RNAseq.g4586.t1       17       98       114       45
                                              101112L9 101112M14 101112M4 101112M9
    Montipora_capitata_HIv3___RNAseq.g4581.t1      360       215      199      386
    Montipora_capitata_HIv3___RNAseq.g4582.t1     1195      1311      431     1076
    Montipora_capitata_HIv3___RNAseq.g4583.t1     2367      1875     1176     2230
    Montipora_capitata_HIv3___RNAseq.g4584.t1      668      1963      800      706
    Montipora_capitata_HIv3___RNAseq.g4585.t1      190       350      103      198
    Montipora_capitata_HIv3___RNAseq.g4586.t1       89       109       18      155
                                              123C14 123C4 123C9 123H14 123H4 123H9
    Montipora_capitata_HIv3___RNAseq.g4581.t1     63   714   351    138   384   146
    Montipora_capitata_HIv3___RNAseq.g4582.t1    643   496  1064    927   378   864
    Montipora_capitata_HIv3___RNAseq.g4583.t1   1349  1591  1529   1377   648  1902
    Montipora_capitata_HIv3___RNAseq.g4584.t1   2724   703   648   3178   497   438
    Montipora_capitata_HIv3___RNAseq.g4585.t1    641    25   146    426     9   208
    Montipora_capitata_HIv3___RNAseq.g4586.t1     78    51   114    121    21    32
                                              123L14 123L4 123L9 123M4 123M9
    Montipora_capitata_HIv3___RNAseq.g4581.t1    190   271   304   771   154
    Montipora_capitata_HIv3___RNAseq.g4582.t1    884   270   838   650   516
    Montipora_capitata_HIv3___RNAseq.g4583.t1   2007   630  1572  1258   892
    Montipora_capitata_HIv3___RNAseq.g4584.t1   1867   330   366   677   397
    Montipora_capitata_HIv3___RNAseq.g4585.t1    368    31    87    15    90
    Montipora_capitata_HIv3___RNAseq.g4586.t1    116     8    62    64    40
                                              131415C14 131415C4 131415C9 131415H14
    Montipora_capitata_HIv3___RNAseq.g4581.t1       198      321      371       196
    Montipora_capitata_HIv3___RNAseq.g4582.t1      1184      712     1338      1314
    Montipora_capitata_HIv3___RNAseq.g4583.t1      2175     1667     2190      1932
    Montipora_capitata_HIv3___RNAseq.g4584.t1      1556      510      473      2267
    Montipora_capitata_HIv3___RNAseq.g4585.t1       310       18      208       453
    Montipora_capitata_HIv3___RNAseq.g4586.t1       143       39      154        82
                                              131415H4 131415H9 131415L14 131415L9
    Montipora_capitata_HIv3___RNAseq.g4581.t1      221      410       170      362
    Montipora_capitata_HIv3___RNAseq.g4582.t1      382     1117      1167     1500
    Montipora_capitata_HIv3___RNAseq.g4583.t1      716     2219      1817     2102
    Montipora_capitata_HIv3___RNAseq.g4584.t1      334      326      2772      534
    Montipora_capitata_HIv3___RNAseq.g4585.t1        0      106       434      129
    Montipora_capitata_HIv3___RNAseq.g4586.t1       15       65       154      144
                                              131415M14 131415M4 131415M9 13M14
    Montipora_capitata_HIv3___RNAseq.g4581.t1       291      271      396    71
    Montipora_capitata_HIv3___RNAseq.g4582.t1      1376      667     1383   912
    Montipora_capitata_HIv3___RNAseq.g4583.t1      2108     1395     2314  1639
    Montipora_capitata_HIv3___RNAseq.g4584.t1      2296      404      492  2194
    Montipora_capitata_HIv3___RNAseq.g4585.t1       375       48      117   527
    Montipora_capitata_HIv3___RNAseq.g4586.t1       156       87       90    92
                                              456C4 456H4 456H9 456L14 456L4 456L9
    Montipora_capitata_HIv3___RNAseq.g4581.t1   196   351   236    125    73   274
    Montipora_capitata_HIv3___RNAseq.g4582.t1   306   458  1301   1306   249  1105
    Montipora_capitata_HIv3___RNAseq.g4583.t1   651   987  1882   1882  1050  2506
    Montipora_capitata_HIv3___RNAseq.g4584.t1   291   547   488   1736   169   495
    Montipora_capitata_HIv3___RNAseq.g4585.t1    29    11    87    384    35   198
    Montipora_capitata_HIv3___RNAseq.g4586.t1    55    29   117    113     0    77
                                              456M14 456M4 456M9 45C14 45C9 45H14
    Montipora_capitata_HIv3___RNAseq.g4581.t1    151   454   250    50  348   411
    Montipora_capitata_HIv3___RNAseq.g4582.t1   1270   525   967   420  981   590
    Montipora_capitata_HIv3___RNAseq.g4583.t1   1649   940  2005   703 2230  1807
    Montipora_capitata_HIv3___RNAseq.g4584.t1   1698   522   595  1402  913  1021
    Montipora_capitata_HIv3___RNAseq.g4585.t1    257    28   154   289  137   238
    Montipora_capitata_HIv3___RNAseq.g4586.t1     72    30    65    17   50    13
                                              67C9 6C14 789H14 789H4 789H9 789L14
    Montipora_capitata_HIv3___RNAseq.g4581.t1  290  114    150   516   332    168
    Montipora_capitata_HIv3___RNAseq.g4582.t1  973 1022    862   628  1338   1027
    Montipora_capitata_HIv3___RNAseq.g4583.t1 1995 1769   1633  1128  2429   1616
    Montipora_capitata_HIv3___RNAseq.g4584.t1  310 1435   2783   963   783   2549
    Montipora_capitata_HIv3___RNAseq.g4585.t1  167  169    489    21   158    457
    Montipora_capitata_HIv3___RNAseq.g4586.t1   77   81     75    68   172     83
                                              789L4 789L9 789M14 789M4 7C14 89C14
    Montipora_capitata_HIv3___RNAseq.g4581.t1   148   409    150   461  208   235
    Montipora_capitata_HIv3___RNAseq.g4582.t1   445  1326    959   551  711  1563
    Montipora_capitata_HIv3___RNAseq.g4583.t1   819  2344   1757  1111 1874  2623
    Montipora_capitata_HIv3___RNAseq.g4584.t1   313   686   2913   906 1307  2195
    Montipora_capitata_HIv3___RNAseq.g4585.t1    34   130    501    30  267   355
    Montipora_capitata_HIv3___RNAseq.g4586.t1    16    82     93    51   55    98
                                              89C9 89M9
    Montipora_capitata_HIv3___RNAseq.g4581.t1  273  230
    Montipora_capitata_HIv3___RNAseq.g4582.t1 1349  823
    Montipora_capitata_HIv3___RNAseq.g4583.t1 2594 1646
    Montipora_capitata_HIv3___RNAseq.g4584.t1  778  624
    Montipora_capitata_HIv3___RNAseq.g4585.t1  271  166
    Montipora_capitata_HIv3___RNAseq.g4586.t1  105   42

- The count matrix shows counts of whole numbers (integers)
- Each column is a sample
- Each row is a gene
- A sample id describes parental crosses (the first string of numbers),
  the pollution exposure (C, L, M, H), and the embryonic stage (4, 9,
  14)

> [!NOTE]
>
> `gcm_filtor` is the gene count matrix that has been filtered down to
> only contain rows where counts are at least 15 in 85% of the samples
> and 2 outlier samples ( `131415L4` and `789C4` ) that had poor HISAT
> alignment and did not group in the exploratory PCA are already
> removed. This step was completed in `code/05_explore`

## Import metadata

``` r
# Load metadata (outlier samples removed, 61 samples remain)
metadata_or <- read_csv(metadata_path)

# set factors
metadata_or <- metadata_or %>% 
  mutate(
    collection_date = as.Date(collection_date, format = "%d-%b-%Y"),
    stage    = factor(stage,    levels = c("cleavage", "prawnchip", "earlygastrula")),
    leachate = factor(leachate, levels = c("control", "low", "mid", "high")),
    spawn_night = factor(spawn_night, levels  = c("July_6th", "July_7th", "July_8th"))
    )

# View metadata structure
str(metadata_or)
```

    tibble [61 × 9] (S3: tbl_df/tbl/data.frame)
     $ sample_id      : chr [1:61] "101112C14" "101112C4" "101112C9" "101112H14" ...
     $ collection_date: Date[1:61], format: "2024-07-08" "2024-07-08" ...
     $ parents        : num [1:61] 101112 101112 101112 101112 101112 ...
     $ group          : chr [1:61] "C14" "C4" "C9" "H14" ...
     $ hpf            : num [1:61] 14 4 9 14 4 9 14 4 9 14 ...
     $ stage          : Factor w/ 3 levels "cleavage","prawnchip",..: 3 1 2 3 1 2 3 1 2 3 ...
     $ leachate       : Factor w/ 4 levels "control","low",..: 1 1 1 4 4 4 2 2 2 3 ...
     $ leachate_mgL   : num [1:61] 0 0 0 1 1 1 0.01 0.01 0.01 0.1 ...
     $ spawn_night    : Factor w/ 3 levels "July_6th","July_7th",..: 3 3 3 3 3 3 3 3 3 3 ...

> [!NOTE]
>
> `metadata_or` contains metadata for 61 samples (*o*utliers *r*emoved)
> with explanatory variables: - `stage` is a factor with 3 levels
> (cleavage, prawnchip, earlygastrula) - `hpf` is a numerical variable
> representing the hours post fertilization 4, 9, and 14 (corresponding
> to `stage`)
>
> - `leachate` is a factor with 4 levels (control, low, mid, high)
>
> - `leachate_mgL` is a continuous numerical variable representing
>   leachate concentrations in mg/L with values 0 , 0.01 , 0.1 , and 1
>   (corresponding to `leachate`)
>
> - `collection_date` is a date-formatted **“%d-%b-%Y”** column with the
>   date of the night embryos were spawned
>
> - `spawn_night` is a factor character string of the spawn night date
>   (either July 6th, July 7th, or July 8th)

# LRT

> [Likelihood ratio
> test](https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#likelihood-ratio-test:~:text=seq%20workflow.-,Likelihood%20ratio%20test,-DESeq2%20offers%20two):
> DESeq2 offers two kinds of hypothesis tests: the Wald test, where we
> use the estimated standard error of a log2 fold change to test if it
> is equal to zero, and the likelihood ratio test (LRT). The LRT
> examines two models for the counts, a full model with a certain number
> of terms and a reduced model, in which some of the terms of the full
> model are removed. The test determines if the increased likelihood of
> the data using the extra terms in the full model is more than expected
> if those extra terms are truly zero. The LRT is therefore useful for
> testing multiple terms at once, for example testing 3 or more levels
> of a factor at once, or all interactions between two variables. The
> LRT for count data is conceptually similar to an analysis of variance
> (ANOVA) calculation in linear regression, except that in the case of
> the Negative Binomial GLM, we use an analysis of deviance (ANODEV),
> where the deviance captures the difference in likelihood between a
> full and a reduced model. The likelihood ratio test can be performed
> by specifying test=“LRT” when using the DESeq function, and providing
> a reduced design formula, e.g. one in which a number of terms from
> design(dds) are removed. The degrees of freedom for the test is
> obtained from the difference between the number of parameters in the
> two models.

> For a control and treatment time series, one can use a design formula
> containing the condition factor, the time factor, and the interaction
> of the two. *In this case, using the likelihood ratio test with a
> reduced model which does not contain the interaction terms will test
> whether the condition induces a change in gene expression at any time
> point after the reference level time point (time 0). An example of the
> later analysis is provided in our [RNA-seq
> workflow](https://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html)*.

We run a few LRT tests :

- full model: `~ spawn_night + stage + leachate + stage:leachate`

- reduced model: `~ spawn_night + stage + leachate` (effect of
  interaction) **INT**

- reduced model: `~ spawn_night + stage` (effects of leachate and stage
  specific leachate together) **JOINT**

&

- full model: `~ spawn_night + stage + leachate`

- reduced model: `~ spawn_night + stage` (tests for main effect of
  leachate) **MAIN**

This tells you if there are any interaction effects between stage and
leachate, but it won’t tell you which interaction is significant, or in
which direction — It helps you decide which model to use. If there
aren’t any differences in significant gene results between models, use
the simpler model.

> [!CAUTION]
>
> Order matters! The first term (`stage`) is what the model controls
> for, and the second term (`leachate`) is the main effect. Meaning the
> model adjusts for stage first and then assesses the leachate effect.

Stage & Leachate will be treated as categorical variables

When embryonic `stage` is treated as a categorical variable DESeq2 will
estimate expression differences between stages, rather than assuming a
linear progression from 4 → 14 hpf. DESeq2 will fit different intercepts
for each stage (i.e., different baseline expression levels). If you
encode stage as a continuous numeric, the model would assume a linear
change in expression per hour of development, which is biologically
risky unless you expect linearity (note: we do not).

> [!NOTE]
>
> `stage` will be treated as a categorical variable

Whether `leachate` is treated as a factor or a numeric (continuous)
variable radically changes the interpretation of both the main effect
and the interaction. If a variable is numeric, DESeq2 treats it as a
continuous covariate. This fits a regression slope, estimating the
change in gene expression per unit increase in the variable. Coding
`leachate` as a categorical factor allows each stage to have its own
independent leachate effect. ✅ Use this if you expect nonlinear,
stage-specific effects of exposure (note: we do).

# Make DESeqDataSet

``` r
dds <- DESeqDataSetFromMatrix(countData = gcm_filtor,
                              colData = metadata_or,
                              design = formula( ~ spawn_night + stage + leachate + stage:leachate))
```

dds \<- DESeqDataSetFromMatrix(countData = gcm_filtor,

colData = metadata_or,

design = formula( ~ spawn_night + stage + leachate + stage:leachate))

## Run DESeq LRT Tests

### Spawn night

Here we treat spawn night as a covariate as a way to control for it as a
nuisance factor.

``` r
dds_sp_LRT <- DESeq(dds, test = "LRT", reduced = ~ stage + leachate + stage:leachate)
```

| Model name     | Reduced model formula                 | What it tests for                         | Interpretation of significant genes                                                                                             |
|----------------|---------------------------------------|-------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------|
| **dds_sp_LRT** | `~ stage + leachate + stage:leachate` | Tests the contribution of **spawn_night** | Genes whose expression significantly varies **among spawn nights**, after accounting for stage, leachate, and their interaction |

``` r
res_spawn <- results(dds_sp_LRT, alpha = 0.05)
summary(res_spawn)
```


    out of 12221 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0 (up)       : 1874, 15%
    LFC < 0 (down)     : 1677, 14%
    outliers [1]       : 6, 0.049%
    low counts [2]     : 0, 0%
    (mean count < 17)
    [1] see 'cooksCutoff' argument of ?results
    [2] see 'independentFiltering' argument of ?results

So this test asks:

“Does including spawn_night explain additional variation in gene
expression, after accounting for stage and leachate (and their
interaction)?”

The fact that you found ~3,551 genes (≈29%) significant means that
spawn_night contributes meaningfully to expression variability. That’s
strong evidence that it is a real batch or parental effect — not random
noise. You’re not interested in biological differences among spawn
nights; you just want to adjust for them so they don’t confound your
leachate comparisons. You include spawn_night in the design formula to
control for it statistically. You do not test it (no LRT against it; no
results table for it). You interpret the coefficients (or LRT results)
for leachate, stage, and stage:leachate — these now represent effects
after adjusting for spawn_night.

### Interaction term

``` r
dds_int_LRT <- DESeq(dds, test = "LRT", reduced = ~ spawn_night + stage + leachate)
```

| Model name      | Reduced model formula              | What it tests for                                              | Interpretation of significant genes                                                               |
|-----------------|------------------------------------|----------------------------------------------------------------|---------------------------------------------------------------------------------------------------|
| **dds_int_LRT** | `~ spawn_night + stage + leachate` | Tests the contribution of the **stage × leachate interaction** | Genes whose expression patterns differ across **leachate treatments in a stage-dependent manner** |

``` r
res_int <- results(dds_int_LRT, alpha = 0.05)
summary(res_int)
```


    out of 12221 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0 (up)       : 6, 0.049%
    LFC < 0 (down)     : 7, 0.057%
    outliers [1]       : 6, 0.049%
    low counts [2]     : 0, 0%
    (mean count < 17)
    [1] see 'cooksCutoff' argument of ?results
    [2] see 'independentFiltering' argument of ?results

Out of 12,221 genes tested:

- Only 13 genes (≈ 0.1 %) were significant (FDR \< 0.05).
- 6 up-regulated, 7 down-regulated.
- Tests **interaction only** (`stage:leachate`).
- You got **13 DEGs** ⇒ This subset show stage-specific responses.

For nearly all genes, the relationship between leachate level and
expression does not change significantly across stages — the leachate
effect is fairly consistent through development. So the interaction term
adds very little explanatory power after controlling for spawn night and
main effects.

This implies `leachate` affects gene expression, but the
direction/magnitude of that effect doesn’t vary much among embryonic
stages (4, 9, 14 hpf).

If you were expecting strong stage-specific responses, this result
suggests that such effects are **minimal** or **limited to a small
subset of genes**.

### Main leachate effect

``` r
dds_main <- DESeqDataSetFromMatrix(countData = gcm_filtor,
                              colData = metadata_or,
                              design = formula( ~ spawn_night + stage + leachate))
```

``` r
dds_main_LRT <- DESeq(dds_main, test = "LRT", reduced = ~ spawn_night + stage)
```

``` r
res_main <- results(dds_main_LRT, alpha = 0.05)
summary(res_main)
```


    out of 12221 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0 (up)       : 0, 0%
    LFC < 0 (down)     : 0, 0%
    outliers [1]       : 14, 0.11%
    low counts [2]     : 0, 0%
    (mean count < 17)
    [1] see 'cooksCutoff' argument of ?results
    [2] see 'independentFiltering' argument of ?results

### Joint

``` r
dds_joint_LRT <- DESeq(dds, test = "LRT", reduced = ~ spawn_night + stage)
```

``` r
res_joint <- results(dds_joint_LRT, alpha = 0.05)
summary(res_joint)
```


    out of 12221 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0 (up)       : 28, 0.23%
    LFC < 0 (down)     : 15, 0.12%
    outliers [1]       : 6, 0.049%
    low counts [2]     : 0, 0%
    (mean count < 17)
    [1] see 'cooksCutoff' argument of ?results
    [2] see 'independentFiltering' argument of ?results

→ Tests the joint contribution of all leachate terms (main effect and
stage:leachate).

→ We got 43 DEGs ⇒ there is leachate-associated signal somewhere (main
and/or interaction taken together).

This indicates some genes (30) have **weak main + weak interaction**
effects that don't cross FDR alone but **do** when evaluated together.

> [!IMPORTANT]
>
> In summary: results are consistent with stage-specific leachate
> responses (interaction), not a uniform, stage-averaged leachate shift
> (main effect). Hence: 0 (main-only), 13 (interaction-only), 43 (any
> leachate signal jointly) is coherent.

> When one performs a likelihood ratio test, the *p* values and the test
> statistic (the `stat` column) are values for the test that removes all
> of the variables which are present in the full design and not in the
> reduced design. This tests the null hypothesis that all the
> coefficients from these variables and levels of these factors are
> equal to zero.
>
> The likelihood ratio test *p* values therefore represent a test of
> *all the variables and all the levels of factors* which are among
> these variables. However, the results table only has space for one
> column of log fold change, so a single variable and a single
> comparison is shown (among the potentially multiple log fold changes
> which were tested in the likelihood ratio test). This is indicated at
> the top of the results table with the text, e.g., log2 fold change
> (MLE): condition C vs A, followed by, LRT p-value: ‘~ batch +
> condition’ vs ‘~ batch’. This indicates that the *p* value is for the
> likelihood ratio test of *all the variables and all the levels*, while
> the log fold change is a single comparison from among those variables
> and levels.

``` r
?results
```

> [Building the results
> table](https://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#multiple-testing:~:text=5.2-,Building%20the%20results%20table,-Calling%20results%20without):
> Calling results without any arguments will extract the estimated log2
> fold changes and p values for the last variable in the design formula.
> If there are more than 2 levels for this variable, results will
> extract the results table for a comparison of the last level over the
> first level.
>
> Genes with small *p* values from the LRT test are those which at one
> or more time points after time 0 showed a strain-specific effect. Note
> therefore that this will not give small *p* values to genes that moved
> up or down over time in the same way in both strains.

``` r
sig_int <- res_int %>%
  as.data.frame() %>%
  filter(padj < 0.05) %>% 
  rownames_to_column(var = "gene_id")

sig_joint <- res_joint%>%
  as.data.frame() %>%
  filter(padj < 0.05) %>% 
  rownames_to_column(var = "gene_id")

overlap_ids <- intersect(sig_int$gene_id, sig_joint$gene_id)
length(overlap_ids)        # how many
```

    [1] 13

``` r
overlap_ids                # which ones
```

     [1] "Montipora_capitata_HIv3___TS.g28558.t3a"    
     [2] "Montipora_capitata_HIv3___TS.g9507.t1"      
     [3] "Montipora_capitata_HIv3___TS.g9699.t2"      
     [4] "Montipora_capitata_HIv3___RNAseq.g40764.t1" 
     [5] "Montipora_capitata_HIv3___TS.g48661.t1"     
     [6] "Montipora_capitata_HIv3___TS.g45793.t1b"    
     [7] "Montipora_capitata_HIv3___RNAseq.g20824.t1a"
     [8] "Montipora_capitata_HIv3___RNAseq.g28195.t1" 
     [9] "Montipora_capitata_HIv3___TS.g15697.t1"     
    [10] "Montipora_capitata_HIv3___RNAseq.g47569.t1" 
    [11] "Montipora_capitata_HIv3___RNAseq.g39262.t1" 
    [12] "Montipora_capitata_HIv3___RNAseq.g18178.t1" 
    [13] "Montipora_capitata_HIv3___RNAseq.g35580.t1" 

> [!NOTE]
>
> The 13 genes found to have stage specific effects (interaction term)
> are also present in the joint list that catches subtle effects of both
> interaction and main leachate. We have 43 genes with significant
> leachate effects overall. This is subtle (43 out of ~12K genes may be
> p-hacked…) but significant based on adjusted pvalues all less than
> 0.05.

> The 43 genes showing significant stage-specific responses to leachate
> exhibit nonlinear or non-monotonic patterns, rather than a simple
> linear trend across leachate levels.

# Wald

> Wald tests assess whether the coefficient **differs significantly from
> 0**

Run the Wald test on the dds object.

``` r
dds <- DESeq(dds, test = "Wald")
```

``` r
resultsNames(dds)
```

     [1] "Intercept"                        "spawn_night_July_7th_vs_July_6th"
     [3] "spawn_night_July_8th_vs_July_6th" "stage_prawnchip_vs_cleavage"     
     [5] "stage_earlygastrula_vs_cleavage"  "leachate_low_vs_control"         
     [7] "leachate_mid_vs_control"          "leachate_high_vs_control"        
     [9] "stageprawnchip.leachatelow"       "stageearlygastrula.leachatelow"  
    [11] "stageprawnchip.leachatemid"       "stageearlygastrula.leachatemid"  
    [13] "stageprawnchip.leachatehigh"      "stageearlygastrula.leachatehigh" 

## Results

**We’re focusing only on leachate effects in this analysis**

### Interaction effects

- “stageprawnchip.leachatelow”

- “stageearlygastrula.leachatelow”

- “stageprawnchip.leachatemid”

- “stageearlygastrula.leachatemid”

- “stageprawnchip.leachatehigh”

- “stageearlygastrula.leachatehigh”

Each of the above interaction terms tells us the extra effect of the
leachate level at a stage (vs the additive sum of the stage and leachate
effects) in comparison to a control baseline count. For the interaction
terms, we are trying to get the “difference of the differences”.

``` r
# prawnchip interaction terms
res_prawnchip_high.int <- results(dds, name = "stageprawnchip.leachatehigh", alpha = 0.05)
res_prawnchip_mid.int <- results(dds, name = "stageprawnchip.leachatemid", alpha = 0.05)
res_prawnchip_low.int  <- results(dds, name = "stageprawnchip.leachatelow", alpha = 0.05)

# early gastrula interaction terms
res_earlygastrula_high.int <- results(dds, name = "stageearlygastrula.leachatehigh", alpha = 0.05)
res_earlygastrula_mid.int <- results(dds, name = "stageearlygastrula.leachatemid", alpha = 0.05)
res_earlygastrula_low.int <- results(dds, name = "stageearlygastrula.leachatelow", alpha = 0.05)
```

#### Prawnchip

##### high

``` r
summary(res_prawnchip_high.int)
```


    out of 12221 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0 (up)       : 0, 0%
    LFC < 0 (down)     : 0, 0%
    outliers [1]       : 6, 0.049%
    low counts [2]     : 0, 0%
    (mean count < 17)
    [1] see 'cooksCutoff' argument of ?results
    [2] see 'independentFiltering' argument of ?results

> [!NOTE]
>
> 0 genes showed a significant stage × leachate interaction in the
> prawnchip stage under the high treatment (Wald test, BH-adjusted p \<
> 0.05).

<div class="column-margin">

> DESeq2 uses size-factor normalization (the “median-ratio method”). It
> corrects for library size and composition so you can compare counts
> across samples. How it works (median-ratio method) Per-gene geometric
> mean across samples (skip genes with all zeros). Per-sample ratios:
> for each gene, divide a sample’s count by that gene’s geometric mean.
> Size factor: the median of those ratios within each sample. Normalized
> counts: normalized = raw_count / size_factor. These size factors are
> used as an offset in the negative binomial GLM for the Wald/LRT tests.
> Normalization is independent of the design formula.

</div>

##### mid

``` r
summary(res_prawnchip_mid.int)
```


    out of 12221 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0 (up)       : 0, 0%
    LFC < 0 (down)     : 0, 0%
    outliers [1]       : 6, 0.049%
    low counts [2]     : 0, 0%
    (mean count < 17)
    [1] see 'cooksCutoff' argument of ?results
    [2] see 'independentFiltering' argument of ?results

> [!NOTE]
>
> 0 genes were found to have significant prawnchip stage x mid leachate
> interaction (adjusted p value \< 0.05)

##### low

``` r
summary(res_prawnchip_low.int)
```


    out of 12221 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0 (up)       : 6, 0.049%
    LFC < 0 (down)     : 1, 0.0082%
    outliers [1]       : 6, 0.049%
    low counts [2]     : 0, 0%
    (mean count < 17)
    [1] see 'cooksCutoff' argument of ?results
    [2] see 'independentFiltering' argument of ?results

> [!NOTE]
>
> 7 genes were found to have significant prawnchip stage x low leachate
> interaction (adjusted p value \< 0.05)

``` r
sig_prawnchip_low.int <- res_prawnchip_low.int %>%
  as.data.frame() %>%
  filter(padj < 0.05) %>% 
  rownames_to_column(var = "gene_id")

write_csv(sig_prawnchip_low.int, file.path(output_path, "sig_prawnchip_low.int.csv"))
```

#### Early gastrula

##### high

``` r
summary(res_earlygastrula_high.int)
```


    out of 12221 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0 (up)       : 0, 0%
    LFC < 0 (down)     : 1, 0.0082%
    outliers [1]       : 6, 0.049%
    low counts [2]     : 0, 0%
    (mean count < 17)
    [1] see 'cooksCutoff' argument of ?results
    [2] see 'independentFiltering' argument of ?results

> [!NOTE]
>
> 1 gene was found to have a significant interaction effect between
> early gastrula and high leachate

``` r
sig_earlygastrula_high.int <- res_earlygastrula_high.int %>%
  as.data.frame() %>%
  filter(padj < 0.05) %>% 
  rownames_to_column(var = "gene_id")

write_csv(sig_earlygastrula_high.int, file.path(output_path, "sig_earlygastrula_high.int.csv"))
```

##### mid

``` r
summary(res_earlygastrula_mid.int)
```


    out of 12221 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0 (up)       : 0, 0%
    LFC < 0 (down)     : 0, 0%
    outliers [1]       : 6, 0.049%
    low counts [2]     : 0, 0%
    (mean count < 17)
    [1] see 'cooksCutoff' argument of ?results
    [2] see 'independentFiltering' argument of ?results

> [!NOTE]
>
> 0 genes were found to have a significant interaction between early
> gastrula and mid leachate

##### low

``` r
summary(res_earlygastrula_low.int)
```


    out of 12221 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0 (up)       : 5, 0.041%
    LFC < 0 (down)     : 3, 0.025%
    outliers [1]       : 6, 0.049%
    low counts [2]     : 0, 0%
    (mean count < 17)
    [1] see 'cooksCutoff' argument of ?results
    [2] see 'independentFiltering' argument of ?results

> [!NOTE]
>
> 8 genes were found to have significant interaction between early
> gastrula and low leachate

``` r
sig_earlygastrula_low.int <- res_earlygastrula_low.int %>%
  as.data.frame() %>%
  filter(padj < 0.05) %>% 
  rownames_to_column(var = "gene_id")

write_csv(sig_earlygastrula_low.int, file.path(output_path, "sig_earlygastrula_low.int.csv"))
```

## Contrasts

See [A guide to designs and contrasts in
DESeq2](https://www.atakanekiz.com/technical/a-guide-to-designs-and-contrasts-in-DESeq2/)
by Dr. Atakan Ekiz

``` r
resultsNames(dds)
```

     [1] "Intercept"                        "spawn_night_July_7th_vs_July_6th"
     [3] "spawn_night_July_8th_vs_July_6th" "stage_prawnchip_vs_cleavage"     
     [5] "stage_earlygastrula_vs_cleavage"  "leachate_low_vs_control"         
     [7] "leachate_mid_vs_control"          "leachate_high_vs_control"        
     [9] "stageprawnchip.leachatelow"       "stageearlygastrula.leachatelow"  
    [11] "stageprawnchip.leachatemid"       "stageearlygastrula.leachatemid"  
    [13] "stageprawnchip.leachatehigh"      "stageearlygastrula.leachatehigh" 

``` r
# cleavage 
res_cleavage_low.vs.control <- results(dds, name = "leachate_low_vs_control", alpha = 0.05)
res_cleavage_mid.vs.control  <- results(dds, name = "leachate_mid_vs_control", alpha = 0.05)
res_cleavage_high.vs.control <- results(dds, name = "leachate_high_vs_control", alpha = 0.05)

# prawnchip
res_prawnchip_low.vs.control <- results(dds, contrast = list(c("leachate_low_vs_control", "stageprawnchip.leachatelow")), alpha = 0.05)
res_prawnchip_mid.vs.control  <- results(dds, contrast = list(c("leachate_mid_vs_control",  "stageprawnchip.leachatemid")), alpha = 0.05)
res_prawnchip_high.vs.control <- results(dds, contrast = list(c("leachate_high_vs_control", "stageprawnchip.leachatehigh")), alpha = 0.05)

# earlygastrula
res_earlygastrula_low.vs.control <- results(dds, contrast = list(c("leachate_low_vs_control", "stageearlygastrula.leachatelow")), alpha = 0.05)
res_earlygastrula_mid.vs.control  <- results(dds, contrast = list(c("leachate_mid_vs_control",  "stageearlygastrula.leachatemid")), alpha = 0.05)
res_earlygastrula_high.vs.control <- results(dds, contrast = list(c("leachate_high_vs_control", "stageearlygastrula.leachatehigh")), alpha = 0.05)
```

### Cleavage

#### high

``` r
summary(res_cleavage_high.vs.control)
```


    out of 12221 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0 (up)       : 2, 0.016%
    LFC < 0 (down)     : 5, 0.041%
    outliers [1]       : 6, 0.049%
    low counts [2]     : 0, 0%
    (mean count < 17)
    [1] see 'cooksCutoff' argument of ?results
    [2] see 'independentFiltering' argument of ?results

> [!NOTE]
>
> 7 genes were found to be significantly differentially expressed
> (adjusted p value \< 0.05) at the cleavage stage in the high leachate
> treatment compared to the control.

``` r
sig_cleavage_high.vs.control <- res_cleavage_high.vs.control %>%
  as.data.frame() %>%
  filter(padj < 0.05) %>% 
  rownames_to_column(var = "gene_id")

write_csv(sig_cleavage_high.vs.control, file.path(output_path, "sig_cleavage_high.vs.control.csv"))
```

#### mid

mid leachate vs. control (mvc) for cleavage

``` r
summary(res_cleavage_mid.vs.control)
```


    out of 12221 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0 (up)       : 57, 0.47%
    LFC < 0 (down)     : 1, 0.0082%
    outliers [1]       : 6, 0.049%
    low counts [2]     : 0, 0%
    (mean count < 17)
    [1] see 'cooksCutoff' argument of ?results
    [2] see 'independentFiltering' argument of ?results

> [!NOTE]
>
> 58 genes were found to be significantly differentially expressed
> (adjusted p value \< 0.05) at the cleavage stage in the mid leachate
> treatment compared to the control. 57 LFC \> 0. 57/58 total genes were
> found in more abundance compared to the controls.

``` r
sig_cleavage_mid.vs.control <- res_cleavage_mid.vs.control %>%
  as.data.frame() %>%
  filter(padj < 0.05) %>% 
  rownames_to_column(var = "gene_id")

write_csv(sig_cleavage_mid.vs.control, file.path(output_path, "sig_cleavage_mid.vs.control.csv"))
```

#### low

for cleavage (baseline stage) \> Because you fit ~ 0 + stage +
leachate + stage:leachate, **the coefficient leachatelow is the low vs
control effect at the reference stage**. From your resultsNames() (no
stagecleavage.leachatelow term is present) the reference stage is
cleavage, so Low vs control at cleavage is simply the leachatelow
coefficient.

``` r
summary(res_cleavage_low.vs.control)
```


    out of 12221 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0 (up)       : 13, 0.11%
    LFC < 0 (down)     : 21, 0.17%
    outliers [1]       : 6, 0.049%
    low counts [2]     : 237, 1.9%
    (mean count < 37)
    [1] see 'cooksCutoff' argument of ?results
    [2] see 'independentFiltering' argument of ?results

> [!NOTE]
>
> 34 genes were found to be significantly differentially expressed
> (adjusted p value \< 0.05) at the cleavage stage in the low leachate
> treatment compared to the control. 13 LFC \> 0 & 21 LFC \< 0.

``` r
sig_cleavage_low.vs.control <- res_cleavage_low.vs.control %>%
  as.data.frame() %>%
  filter(padj < 0.05) %>% 
  rownames_to_column(var = "gene_id")

write_csv(sig_cleavage_low.vs.control, file.path(output_path, "sig_cleavage_low.vs.control.csv"))
```

### Prawnchip

#### high

for prawnchip

``` r
summary(res_prawnchip_high.vs.control)
```


    out of 12221 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0 (up)       : 0, 0%
    LFC < 0 (down)     : 0, 0%
    outliers [1]       : 6, 0.049%
    low counts [2]     : 0, 0%
    (mean count < 17)
    [1] see 'cooksCutoff' argument of ?results
    [2] see 'independentFiltering' argument of ?results

> [!NOTE]
>
> 0 genes differentially expressed (adjusted p value \< 0.05) at the
> prawn chip stage in the high leachate treatment when compared to the
> control.

#### mid

``` r
summary(res_prawnchip_mid.vs.control)
```


    out of 12221 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0 (up)       : 0, 0%
    LFC < 0 (down)     : 0, 0%
    outliers [1]       : 6, 0.049%
    low counts [2]     : 0, 0%
    (mean count < 17)
    [1] see 'cooksCutoff' argument of ?results
    [2] see 'independentFiltering' argument of ?results

> [!NOTE]
>
> 0 genes differentially expressed (adjusted p value \< 0.05) at the
> prawn chip stage in the high leachate treatment when compared to the
> control.

#### low

``` r
summary(res_prawnchip_low.vs.control)
```


    out of 12221 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0 (up)       : 0, 0%
    LFC < 0 (down)     : 0, 0%
    outliers [1]       : 6, 0.049%
    low counts [2]     : 0, 0%
    (mean count < 17)
    [1] see 'cooksCutoff' argument of ?results
    [2] see 'independentFiltering' argument of ?results

> [!NOTE]
>
> 0 genes differentially expressed (adjusted p value \< 0.05) at the
> prawn chip stage in the high leachate treatment when compared to the
> control.

### Early gastrula

#### high

``` r
summary(res_earlygastrula_high.vs.control)
```


    out of 12221 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0 (up)       : 27, 0.22%
    LFC < 0 (down)     : 0, 0%
    outliers [1]       : 6, 0.049%
    low counts [2]     : 9709, 79%
    (mean count < 1467)
    [1] see 'cooksCutoff' argument of ?results
    [2] see 'independentFiltering' argument of ?results

> [!NOTE]
>
> 27 genes differentially expressed (adjusted p value \< 0.05) at the
> early gastrula stage in the high leachate treatment when compared to
> the control. All 27 of these genes had more abundant transcripts in
> treated samples compared to controls.

``` r
sig_earlygastrula_high.vs.control <- res_earlygastrula_high.vs.control %>%
  as.data.frame() %>%
  filter(padj < 0.05) %>% 
  rownames_to_column(var = "gene_id")

write_csv(sig_earlygastrula_high.vs.control, file.path(output_path, "sig_earlygastrula_high.vs.control.csv"))
```

#### mid

``` r
summary(res_earlygastrula_mid.vs.control)
```


    out of 12221 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0 (up)       : 0, 0%
    LFC < 0 (down)     : 0, 0%
    outliers [1]       : 6, 0.049%
    low counts [2]     : 0, 0%
    (mean count < 17)
    [1] see 'cooksCutoff' argument of ?results
    [2] see 'independentFiltering' argument of ?results

> [!NOTE]
>
> 0 genes differentially expressed (adjusted p value \< 0.05) at the
> early gastrula stage in the mid leachate treatment when compared to
> the control.

#### low

``` r
summary(res_earlygastrula_low.vs.control)
```


    out of 12221 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0 (up)       : 0, 0%
    LFC < 0 (down)     : 0, 0%
    outliers [1]       : 6, 0.049%
    low counts [2]     : 0, 0%
    (mean count < 17)
    [1] see 'cooksCutoff' argument of ?results
    [2] see 'independentFiltering' argument of ?results

> [!NOTE]
>
> 0 genes differentially expressed (adjusted p value \< 0.05) at the
> early gastrula stage in the low leachate treatment when compared to
> the control.

# Summary

So we have 9 different ‘significant gene’ lists. Each list indicates
genes in which leachate had a significant effect.

LRT: - `sig_int`, 13 DEGs - `sig_joint`. 44 DEGs

Interaction effects - `sig_prawnchip_low.int`, 7 DEGs -
`sig_earlygastrula_high.int`, 1 DEGs - `sig_earlygastrula_low.int`, 8
DEGs

Contrasts - `sig_cleavage_low.vs.control`, 34 DEGs -
`sig_cleavage_mid.vs.control`, 58 DEGs - `sig_cleavage_high.vs.control`,
7 DEGs - `sig_earlygastrula_high.vs.control`, 27 DEGs

``` r
# 1) Put all of your sig data frames into a *named* list
sig_lists <- list(
  int                           = sig_int,
  joint                         = sig_joint,
  prawnchip_low.int             = sig_prawnchip_low.int,
  earlygastrula_high.int        = sig_earlygastrula_high.int,
  earlygastrula_low.int         = sig_earlygastrula_low.int,
  cleavage_low.vs.control       = sig_cleavage_low.vs.control,
  cleavage_mid.vs.control       = sig_cleavage_mid.vs.control,
  cleavage_high.vs.control      = sig_cleavage_high.vs.control,
  earlygastrula_high.vs.control = sig_earlygastrula_high.vs.control
)

# 2) Smoosh: add a 'list' column indicating the source and stack them
sig_summary <- bind_rows(sig_lists, .id = "list") %>%
  relocate(list, gene_id)

write_csv(sig_summary, file.path(output_path, "sig_summary.csv"))
```

``` r
sig_all <- dplyr::distinct(sig_summary, gene_id)
nrow(sig_all)
```

    [1] 130

``` r
write_csv(sig_all, file.path(output_path, "sig_all.csv"))
```

There are 130 total DEGs across LRT and Wald tests where leachate popped
up as a significant covariate. \# Save count matrices

Make variance stabilized count matrix of these 130 genes of interest

``` r
# VST transform
vsd <- vst(dds, blind = FALSE) # DESeqTransform Object
vst_mat <- assay(vsd) # VST count matrix with genes as rows and samples as columns

# Keep the 130 genes of interest
vst_mat_130 <- vst_mat[sig_all$gene_id, , drop = FALSE]
dim(vst_mat_130)
```

    [1] 130  61

``` r
# Convert to dataframes and save 
vst_mat_df <- vst_mat %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene_id")

vst_mat_130_df <- vst_mat_130 %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene_id")

write_csv(vst_mat_df, file.path(output_path,"vst_mat_df.csv"))
write_csv(vst_mat_130_df, file.path(output_path,"vst_mat_130_df.csv"))
```

Make library-size adjusted count matrix of these 130 genes of interest

``` r
norm <- counts(dds, normalized = TRUE) # size factor normalized count matrix with genes as rows and samples as columns

# Keep the 130 genes of interest
norm_130 <- norm[sig_all$gene_id, , drop = FALSE]
dim(norm_130)
```

    [1] 130  61

``` r
# Convert to dataframes and save 
norm_df <- norm %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene_id")

norm_130_df <- norm_130 %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene_id")

write_csv(norm_df, file.path(output_path,"norm_df.csv"))
write_csv(norm_130_df, file.path(output_path,"norm_130_df.csv"))
```

Make raw count matrix of these 130 genes of interest

``` r
raw <- counts(dds) # raw count matrix with genes as rows and samples as columns

# Keep the 130 genes of interest
raw_130 <- raw[sig_all$gene_id, , drop = FALSE]
dim(raw_130)
```

    [1] 130  61

``` r
# Convert to dataframes and save 
raw_df <- raw %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene_id")

raw_130_df <- raw_130 %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene_id")

write_csv(raw_df, file.path(output_path,"raw_df.csv"))
write_csv(raw_130_df, file.path(output_path,"raw_130_df.csv"))
```

# Next steps

- Make a PCA and PERMANOVA of the DEGs …
- I’d simply like an easy absence/presence heatmap that shows all the
  genes, and which lists they are significant in.
- I’d then like to show on the heatmap the log2FC, and add astericks for
  significance levels p\<0.05, p\<0.01, p\<0.001.
- I’d like gene count plots visually summarizing each list type (LRT,
  Interaction, Stagexleachate effects) in tiny tidy multiples
- Do these genes cluster into expression patterns? … Clustering
- What do these genes do ? … GO enrichment
- Are there correlations between these differentially expressed genes
  and differentially abundant taxa?
