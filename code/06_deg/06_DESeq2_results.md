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
    - [<span class="toc-section-number">2.4.1</span> Reorder metadata
      factors](#reorder-metadata-factors)
- [<span class="toc-section-number">3</span> LRT](#lrt)
- [<span class="toc-section-number">4</span> Categorical
  vs. continuous](#categorical-vs-continuous)
- [<span class="toc-section-number">5</span> Make
  DESeqDataSet](#make-deseqdataset)
  - [<span class="toc-section-number">5.1</span> Run DESeq LRT
    Tests](#run-deseq-lrt-tests)
  - [<span class="toc-section-number">5.2</span> Results](#results)
    - [<span class="toc-section-number">5.2.1</span>
      categorical](#categorical)
    - [<span class="toc-section-number">5.2.2</span>
      continuous](#continuous)
    - [<span class="toc-section-number">5.2.3</span> Save significant
      LRT genes](#save-significant-lrt-genes)
    - [<span class="toc-section-number">5.2.4</span> Save LRT test
      normalized count matrix](#save-lrt-test-normalized-count-matrix)
- [<span class="toc-section-number">6</span> Wald](#wald)
  - [<span class="toc-section-number">6.0.1</span> Save Wald test
    normalized count matrix](#save-wald-test-normalized-count-matrix)
  - [<span class="toc-section-number">6.1</span> Results](#results-1)
    - [<span class="toc-section-number">6.1.1</span> Interaction
      effects](#interaction-effects)
  - [<span class="toc-section-number">6.2</span> Stage x leachate
    effects](#stage-x-leachate-effects)
    - [<span class="toc-section-number">6.2.1</span>
      Cleavage](#cleavage)
    - [<span class="toc-section-number">6.2.2</span>
      Prawnchip](#prawnchip-1)
    - [<span class="toc-section-number">6.2.3</span> Early
      gastrula](#early-gastrula-1)
- [<span class="toc-section-number">7</span> Summary](#summary)
- [<span class="toc-section-number">8</span> Next steps](#next-steps)

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

- `gcm_filtor.csv` is the filtered gene count matrix for 61 samples (two
  outliers removed).
- `metadata_or.csv` is the metadata for 61 samples (two outliers
  removed)

## Outputs

- `output/06_deg/resul/int_wald_heatmap.png`
- `volcano_facet_interaction.png`

# Setup

``` r
output <- "../../output/06_deg/results/"
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
> only contain rows where counts are at least 15 in 75% of the samples
> and 2 outlier samples ( `131415L4` and `789C4` ) that had poor HISAT
> alignment and did not group in the exploratory PCA are already
> removed. This step was completed in `code/05_explore`

## Import metadata

``` r
metadata_or <- read_csv("../../metadata/metadata_or.csv")
str(metadata_or)
```

    spc_tbl_ [61 × 7] (S3: spec_tbl_df/tbl_df/tbl/data.frame)
     $ sample_id   : chr [1:61] "101112C14" "101112C4" "101112C9" "101112H14" ...
     $ parents     : num [1:61] 101112 101112 101112 101112 101112 ...
     $ group       : chr [1:61] "C14" "C4" "C9" "H14" ...
     $ hpf         : num [1:61] 14 4 9 14 4 9 14 4 9 14 ...
     $ stage       : chr [1:61] "earlygastrula" "cleavage" "prawnchip" "earlygastrula" ...
     $ leachate    : chr [1:61] "control" "control" "control" "high" ...
     $ leachate_mgL: num [1:61] 0 0 0 1 1 1 0.01 0.01 0.01 0.1 ...
     - attr(*, "spec")=
      .. cols(
      ..   sample_id = col_character(),
      ..   parents = col_double(),
      ..   group = col_character(),
      ..   hpf = col_double(),
      ..   stage = col_character(),
      ..   leachate = col_character(),
      ..   leachate_mgL = col_double()
      .. )
     - attr(*, "problems")=<externalptr> 

### Reorder metadata factors

``` r
# Reorder factors
metadata_or$stage <- as.factor(metadata_or$stage)
metadata_or$stage <- fct_relevel(metadata_or$stage, "cleavage", "prawnchip", "earlygastrula")
 
metadata_or$leachate <- as.factor(metadata_or$leachate)
metadata_or$leachate <- fct_relevel(metadata_or$leachate, "control", "low", "mid", "high")

str(metadata_or)
```

    spc_tbl_ [61 × 7] (S3: spec_tbl_df/tbl_df/tbl/data.frame)
     $ sample_id   : chr [1:61] "101112C14" "101112C4" "101112C9" "101112H14" ...
     $ parents     : num [1:61] 101112 101112 101112 101112 101112 ...
     $ group       : chr [1:61] "C14" "C4" "C9" "H14" ...
     $ hpf         : num [1:61] 14 4 9 14 4 9 14 4 9 14 ...
     $ stage       : Factor w/ 3 levels "cleavage","prawnchip",..: 3 1 2 3 1 2 3 1 2 3 ...
     $ leachate    : Factor w/ 4 levels "control","low",..: 1 1 1 4 4 4 2 2 2 3 ...
     $ leachate_mgL: num [1:61] 0 0 0 1 1 1 0.01 0.01 0.01 0.1 ...
     - attr(*, "spec")=
      .. cols(
      ..   sample_id = col_character(),
      ..   parents = col_double(),
      ..   group = col_character(),
      ..   hpf = col_double(),
      ..   stage = col_character(),
      ..   leachate = col_character(),
      ..   leachate_mgL = col_double()
      .. )
     - attr(*, "problems")=<externalptr> 

> [!NOTE]
>
> `metadata_or` contains metadata for 61 samples (*o*utliers *r*emoved)
> with explanatory variables: - `stage` is a factor with 3 levels
> (cleavage, prawnchip, earlygastrula) - `hpf` is a discrete numerical
> variable representing the hours post fertilization 4, 9, and 14
> (corresponding to `stage`) - `leachate` is a factor with 4 levels
> (control, low, mid, high) - `leachate_mgL` is a continuous numerical
> variable representing leachate concentrations in mg/L with values 0 ,
> 0.01 , 0.1 , and 1 (corresponding to `leachate`)

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

This compares:

- full model: `~ stage + leachate + stage:leachate`

- reduced model: `~ stage + leachate`

This tells you if there are any interaction effects between stage and
leachate, but it won’t tell you which interaction is significant, or in
which direction — you need Wald contrasts for that.

> [!CAUTION]
>
> Order matters! The first term (`stage`) is what the model controls
> for, and the second term (`leachate`) is the main effect. Meaning the
> model adjusts for stage first and then assesses the leachate effect.

# Categorical vs. continuous

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
dds_cat <- DESeqDataSetFromMatrix(countData = gcm_filtor,
                              colData = metadata_or,
                              design = formula( ~0 + stage + leachate + stage:leachate))
```

Coding `leachate_mgL` as a continuous numeric tests for a change in
expression per unit increase in leachate, and the interaction term tests
whether the slope of the dose-response differs depending on stage. ✅
Use this if you expect a linear or dose-dependent trend, and want to ask
if the strength or direction of the trend varies by stage (exploring if
this is the case).

``` r
dds_con <- DESeqDataSetFromMatrix(countData = gcm_filtor,
                              colData = metadata_or,
                              design = formula( ~0 + stage + leachate_mgL + stage:leachate_mgL))
```

## Run DESeq LRT Tests

> Does modeling leachate as a factor vs. a continuous dose better
> explain gene expression changes across stages?

> [!NOTE]
>
> ~0 zero-intercept model

``` r
dds_cat_LRT <- DESeq(dds_cat, test = "LRT", reduced = ~0 + stage + leachate)
dds_con_LRT <- DESeq(dds_con, test = "LRT", reduced = ~0 + stage + leachate_mgL)
```

| Model         | Full model                                       | Reduced model               | LRT tests…                                                 |
|---------------|--------------------------------------------------|-----------------------------|------------------------------------------------------------|
| `dds_cat_LRT` | `~0 + stage + leachate + stage:leachate`         | `~0 + stage + leachate`     | Are there **stage-specific effects of leachate (factor)**? |
| `dds_con_LRT` | `~0 + stage + leachate_mgL + stage:leachate_mgL` | `~0 + stage + leachate_mgL` | Are there **stage-specific trends** in dose response?      |

Plot dispersion estimates

``` r
plotDispEsts(dds_cat_LRT)
```

<img
src="06_DESeq2_results_files/figure-commonmark/unnamed-chunk-10-1.jpeg"
width="4200" />

``` r
plotDispEsts(dds_con_LRT)
```

<img
src="06_DESeq2_results_files/figure-commonmark/unnamed-chunk-10-2.jpeg"
width="4200" />

> [!NOTE]
>
> What does the categorical models dispersion (.01 - 1) vs the
> continuous model (.005 - 0.5) show?! Pay attention to the scale of the
> Y axis when comparing them.

## Results

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
> and levels. See the help page for *results* for more details.

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

### categorical

``` r
resultsNames(dds_cat_LRT)
```

     [1] "stagecleavage"                   "stageprawnchip"                 
     [3] "stageearlygastrula"              "leachatelow"                    
     [5] "leachatemid"                     "leachatehigh"                   
     [7] "stageprawnchip.leachatelow"      "stageearlygastrula.leachatelow" 
     [9] "stageprawnchip.leachatemid"      "stageearlygastrula.leachatemid" 
    [11] "stageprawnchip.leachatehigh"     "stageearlygastrula.leachatehigh"

``` r
res_cat_LRT <- results(dds_cat_LRT, alpha = 0.05)
res_cat_LRT
```

    log2 fold change (MLE): stageearlygastrula.leachatehigh 
    LRT p-value: '~ 0 + stage + leachate + stage:leachate' vs '~ 0 + stage + leachate' 
    DataFrame with 12565 rows and 6 columns
                                                baseMean log2FoldChange     lfcSE
                                               <numeric>      <numeric> <numeric>
    Montipora_capitata_HIv3___RNAseq.g4581.t1    297.951      0.8443370  0.465080
    Montipora_capitata_HIv3___RNAseq.g4582.t1    801.344     -0.3244420  0.245480
    Montipora_capitata_HIv3___RNAseq.g4583.t1   1571.805      0.0490377  0.294234
    Montipora_capitata_HIv3___RNAseq.g4584.t1    989.702     -0.3842521  0.371563
    Montipora_capitata_HIv3___RNAseq.g4585.t1    168.821      0.3265516  0.631096
    ...                                              ...            ...       ...
    Montipora_capitata_HIv3___RNAseq.g4570.t1   159.2361      0.2601221  0.497611
    Montipora_capitata_HIv3___RNAseq.g14.t1      60.9601      0.0843803  0.516090
    Montipora_capitata_HIv3___RNAseq.g3903.t1    46.7101     -0.1664749  0.658533
    Montipora_capitata_HIv3___RNAseq.g29883.t1  240.1779     -0.8602191  0.610604
    Montipora_capitata_HIv3___RNAseq.g45863.t1  108.5158     -0.3103637  0.492484
                                                    stat    pvalue      padj
                                               <numeric> <numeric> <numeric>
    Montipora_capitata_HIv3___RNAseq.g4581.t1    7.95135 0.2416894  0.999984
    Montipora_capitata_HIv3___RNAseq.g4582.t1    6.25767 0.3949539  0.999984
    Montipora_capitata_HIv3___RNAseq.g4583.t1    8.02957 0.2359452  0.999984
    Montipora_capitata_HIv3___RNAseq.g4584.t1    2.70745 0.8445671  0.999984
    Montipora_capitata_HIv3___RNAseq.g4585.t1   11.17537 0.0831053  0.999984
    ...                                              ...       ...       ...
    Montipora_capitata_HIv3___RNAseq.g4570.t1    1.74820  0.941340  0.999984
    Montipora_capitata_HIv3___RNAseq.g14.t1      3.76475  0.708476  0.999984
    Montipora_capitata_HIv3___RNAseq.g3903.t1    5.52177  0.478830  0.999984
    Montipora_capitata_HIv3___RNAseq.g29883.t1   4.69604  0.583346  0.999984
    Montipora_capitata_HIv3___RNAseq.g45863.t1   4.89254  0.557668  0.999984

``` r
mcols(res_cat_LRT, use.names = TRUE)
```

    DataFrame with 6 rows and 2 columns
                           type            description
                    <character>            <character>
    baseMean       intermediate mean of normalized c..
    log2FoldChange      results log2 fold change (ML..
    lfcSE               results standard error: stag..
    stat                results LRT statistic: '~ 0 ..
    pvalue              results LRT p-value: '~ 0 + ..
    padj                results   BH adjusted p-values

``` r
summary(res_cat_LRT)
```


    out of 12565 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0 (up)       : 28, 0.22%
    LFC < 0 (down)     : 16, 0.13%
    outliers [1]       : 0, 0%
    low counts [2]     : 0, 0%
    (mean count < 23)
    [1] see 'cooksCutoff' argument of ?results
    [2] see 'independentFiltering' argument of ?results

### continuous

``` r
resultsNames(dds_con_LRT)
```

    [1] "stagecleavage"                   "stageprawnchip"                 
    [3] "stageearlygastrula"              "leachate_mgL"                   
    [5] "stageprawnchip.leachate_mgL"     "stageearlygastrula.leachate_mgL"

``` r
res_con_LRT <- results(dds_con_LRT, alpha = 0.05)
summary(res_con_LRT)
```


    out of 12565 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0 (up)       : 2, 0.016%
    LFC < 0 (down)     : 0, 0%
    outliers [1]       : 1, 0.008%
    low counts [2]     : 0, 0%
    (mean count < 23)
    [1] see 'cooksCutoff' argument of ?results
    [2] see 'independentFiltering' argument of ?results

### Save significant LRT genes

``` r
sig_cat_LRT <- res_cat_LRT %>%
  as.data.frame() %>%
  filter(padj < 0.05) %>% 
  rownames_to_column(var = "gene_id")

sig_con_LRT <- res_con_LRT %>%
  as.data.frame() %>%
  filter(padj < 0.05) %>% 
  rownames_to_column(var = "gene_id")

overlap_ids <- intersect(sig_cat_LRT$gene_id, sig_con_LRT$gene_id)
length(overlap_ids)        # how many
```

    [1] 2

``` r
overlap_ids                # which ones
```

    [1] "Montipora_capitata_HIv3___RNAseq.g11129.t1"
    [2] "Montipora_capitata_HIv3___RNAseq.g21373.t1"

> [!NOTE]
>
> The 2 significant genes (“Montipora_capitata_HIv3\_**RNAseq.g11129.t1”
> & ”Montipora_capitata_HIv3**\_RNAseq.g21373.t1”) found using leachate
> as a continuous factor are also included in the 44 significant genes
> found when leachate was coded as a categorical factor.
>
> This means that 2 out of 44 significant genes may have a linear
> relationship (compared to a non-linear one!) as we look at the
> interactions between leachate and stage. **All 44 of these genes
> respond to leachate differently depending on stage. Response to
> leachate pollution differs by embryonic stage.**
>
> 44 genes show significant stage-by-leachate interaction when leachate
> is treated categorically → meaning the response to leachate varies by
> stage, but not necessarily in a linear/dose-dependent way.
>
> 2 genes show significant stage-specific dose-response effects when
> leachate is treated as continuous → they show a linear or monotonic
> trend across leachate concentrations that varies by stage.
>
> The 2 overlapping genes are well-modeled by a linear trend in leachate
> (i.e., a slope, not just group differences). The 42 other genes show
> non-linear or non-monotonic behavior — meaning their response does not
> follow a consistent increase/decrease with dose. Instead, they respond
> differently to each specific leachate level depending on the stage
> (e.g., they might increase at 0.1 mg/L but decrease at 1 mg/L, or vice
> versa).

Because the 2 genes that are captured in the continuous model are also
caputured in the categorical model we will just save the significant
categorical list

``` r
write_csv(sig_cat_LRT, file.path(output,"sig_cat_LRT.csv"))
```

> [!NOTE]
>
> `sig_cat_LRT` is the DESeq LRT result that shows 44 genes with
> significant interactions between stage and leachate!

> Of the 44 genes showing significant stage-specific responses to
> leachate (based on the categorical model), only 2 also show a
> significant linear dose-response interaction with stage (based on the
> continuous model). This suggests that while stage modulates leachate
> response broadly, most genes (42/44) exhibit nonlinear or
> non-monotonic patterns, rather than a simple linear trend across
> leachate levels.

### Save LRT test normalized count matrix

``` r
# 1. extract normalized counts (matrix)
norm_count_mat_LRT <- counts(dds_cat_LRT, normalized = TRUE)
# column ids are sample_id from metadata_or
# rownames are gene_id
```

> [!NOTE]
>
> We have 44 genes with significant interaction effects overall. Next
> step: 👉 Use Wald tests to dig into which interactions (e.g., stage ×
> leachate level) are responsible — gene by gene.

# Wald

> Wald tests on interaction terms assess whether the interaction
> coefficient **differs significantly from 0** — i.e., whether the
> combined effect of that `hpf` and `leachate` level deviates from what
> would be expected by the additive effects alone.

| What You Want to Know                                               | Use            | Explanation                                                                   |
|---------------------------------------------------------------------|----------------|-------------------------------------------------------------------------------|
| Specific interaction term significant?                              | Wald           | `results(dds, name = "stageX.leachateY")`                                     |
| Total effect of a leachate treatment at a specific stage vs control | Wald contrasts | `results(dds, contrast = list(c("leachateY_vs_control", "stageX.leachateY"))` |

Run the Wald test on the categorical dds object.

``` r
dds_wald <- DESeq(dds_cat, test = "Wald")
```

### Save Wald test normalized count matrix

``` r
# 1. extract normalized counts (matrix)
norm_count_mat_Wald <- counts(dds_wald, normalized = TRUE)
# column ids are sample_id from metadata_or
# rownames are gene_id
```

Are the underlying normalized counts for the Wald and LRT the same?

``` r
identical(norm_count_mat_LRT, norm_count_mat_Wald)
```

    [1] TRUE

Ok so just save one normalized count matrix

``` r
# 1. extract normalized counts (matrix)
norm_count_mat <- counts(dds_wald, normalized = TRUE) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene_id")
# column ids are sample_id from metadata_or
# rownames are gene_id

write_csv(norm_count_mat, file.path(output,"norm_count_mat.csv"))
```

``` r
resultsNames(dds_wald)
```

     [1] "stagecleavage"                   "stageprawnchip"                 
     [3] "stageearlygastrula"              "leachatelow"                    
     [5] "leachatemid"                     "leachatehigh"                   
     [7] "stageprawnchip.leachatelow"      "stageearlygastrula.leachatelow" 
     [9] "stageprawnchip.leachatemid"      "stageearlygastrula.leachatemid" 
    [11] "stageprawnchip.leachatehigh"     "stageearlygastrula.leachatehigh"

## Results

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
terms, we are trying to get the “difference of the differences”. For
this analysis, we are interested in genes that are differentially
expressed in treatment conditions, but not in the control condition (ie.
genes uniquely changing in a stage due to treatment).

``` r
# prawnchip interaction terms
res_prawnchip_high.int <- results(dds_wald, name = "stageprawnchip.leachatehigh", alpha = 0.05)
res_prawnchip_mid.int <- results(dds_wald, name = "stageprawnchip.leachatemid", alpha = 0.05)
res_prawnchip_low.int  <- results(dds_wald, name = "stageprawnchip.leachatelow", alpha = 0.05)

# early gastrula interaction terms
res_earlygastrula_high.int <- results(dds_wald, name = "stageearlygastrula.leachatehigh", alpha = 0.05)
res_earlygastrula_mid.int <- results(dds_wald, name = "stageearlygastrula.leachatemid", alpha = 0.05)
res_earlygastrula_low.int <- results(dds_wald, name = "stageearlygastrula.leachatelow", alpha = 0.05)
```

#### Prawnchip

##### high

``` r
summary(res_prawnchip_high.int)
```


    out of 12565 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0 (up)       : 1, 0.008%
    LFC < 0 (down)     : 0, 0%
    outliers [1]       : 0, 0%
    low counts [2]     : 0, 0%
    (mean count < 23)
    [1] see 'cooksCutoff' argument of ?results
    [2] see 'independentFiltering' argument of ?results

> [!NOTE]
>
> 1 gene was found to have significant stage x leachate interaction
> (adjusted p value \< 0.05)

``` r
sig_prawnchip_high.int <- res_prawnchip_high.int %>%
  as.data.frame() %>%
  filter(padj < 0.05) %>% 
  rownames_to_column(var = "gene_id")

sig_prawnchip_high.int
```

                                         gene_id baseMean log2FoldChange    lfcSE
    1 Montipora_capitata_HIv3___RNAseq.g11129.t1 31.18656       5.762848 1.214606
          stat       pvalue       padj
    1 4.744625 2.088931e-06 0.02624741

> [!NOTE]
>
> Only one gene showed a significant stage × leachate interaction in the
> prawnchip stage under the high treatment (Wald test, BH-adjusted p \<
> 0.05). Specifically,
> \*\*Montipora_capitata_HIv3\_\_\_RNAseq.g11129.t1\*\* was up-regulated
> during the prawnchip stage in the high treatment (log2FC = +5.76, padj
> = 0.026) relative to the expected count levels…

Save significant gene list

``` r
write_csv(sig_prawnchip_high.int, file.path(output, "sig_prawnchip_high.int.csv"))
```

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


    out of 12565 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0 (up)       : 0, 0%
    LFC < 0 (down)     : 5, 0.04%
    outliers [1]       : 0, 0%
    low counts [2]     : 0, 0%
    (mean count < 23)
    [1] see 'cooksCutoff' argument of ?results
    [2] see 'independentFiltering' argument of ?results

> [!NOTE]
>
> 5 genes were found to have significant prawnchip stage x mid leachate
> interaction (adjusted p value \< 0.05)

``` r
sig_prawnchip_mid.int <- res_prawnchip_mid.int %>%
  as.data.frame() %>%
  filter(padj < 0.05) %>% 
  rownames_to_column(var = "gene_id")

write_csv(sig_prawnchip_mid.int, file.path(output, "sig_prawnchip_mid.int.csv"))
```

##### low

``` r
summary(res_prawnchip_low.int)
```


    out of 12565 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0 (up)       : 10, 0.08%
    LFC < 0 (down)     : 2, 0.016%
    outliers [1]       : 0, 0%
    low counts [2]     : 0, 0%
    (mean count < 23)
    [1] see 'cooksCutoff' argument of ?results
    [2] see 'independentFiltering' argument of ?results

> [!NOTE]
>
> 12 genes were found to have significant prawnchip stage x low leachate
> interaction (adjusted p value \< 0.05)

``` r
sig_prawnchip_low.int <- res_prawnchip_low.int %>%
  as.data.frame() %>%
  filter(padj < 0.05) %>% 
  rownames_to_column(var = "gene_id")

write_csv(sig_prawnchip_low.int, file.path(output, "sig_prawnchip_low.int.csv"))
```

#### Early gastrula

##### high

``` r
summary(res_earlygastrula_high.int)
```


    out of 12565 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0 (up)       : 2, 0.016%
    LFC < 0 (down)     : 1, 0.008%
    outliers [1]       : 0, 0%
    low counts [2]     : 0, 0%
    (mean count < 23)
    [1] see 'cooksCutoff' argument of ?results
    [2] see 'independentFiltering' argument of ?results

> [!NOTE]
>
> 3 genes were found to have a significant interaction effect between
> early gastrula and high leachate

``` r
sig_earlygastrula_high.int <- res_earlygastrula_high.int %>%
  as.data.frame() %>%
  filter(padj < 0.05) %>% 
  rownames_to_column(var = "gene_id")

write_csv(sig_earlygastrula_high.int, file.path(output, "sig_earlygastrula_high.int.csv"))
```

##### mid

``` r
summary(res_earlygastrula_mid.int)
```


    out of 12565 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0 (up)       : 0, 0%
    LFC < 0 (down)     : 11, 0.088%
    outliers [1]       : 0, 0%
    low counts [2]     : 5603, 45%
    (mean count < 336)
    [1] see 'cooksCutoff' argument of ?results
    [2] see 'independentFiltering' argument of ?results

> [!NOTE]
>
> 11 genes were found to have a significant interaction between early
> gastrula and mid leachate

``` r
sig_earlygastrula_mid.int <- res_earlygastrula_mid.int %>%
  as.data.frame() %>%
  filter(padj < 0.05) %>% 
  rownames_to_column(var = "gene_id")

write_csv(sig_earlygastrula_mid.int, file.path(output, "sig_earlygastrula_mid.int.csv"))
```

##### low

``` r
summary(res_earlygastrula_low.int)
```


    out of 12565 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0 (up)       : 7, 0.056%
    LFC < 0 (down)     : 4, 0.032%
    outliers [1]       : 0, 0%
    low counts [2]     : 0, 0%
    (mean count < 23)
    [1] see 'cooksCutoff' argument of ?results
    [2] see 'independentFiltering' argument of ?results

> [!NOTE]
>
> 11 genes were found to have significant interaction between early
> gastrula and low leachate

``` r
sig_earlygastrula_low.int <- res_earlygastrula_low.int %>%
  as.data.frame() %>%
  filter(padj < 0.05) %>% 
  rownames_to_column(var = "gene_id")

write_csv(sig_earlygastrula_low.int, file.path(output, "sig_earlygastrula_low.int.csv"))
```

## Stage x leachate effects

See [A guide to designs and contrasts in
DESeq2](https://www.atakanekiz.com/technical/a-guide-to-designs-and-contrasts-in-DESeq2/)
by Dr. Atakan Ekiz

``` r
resultsNames(dds_wald)
```

     [1] "stagecleavage"                   "stageprawnchip"                 
     [3] "stageearlygastrula"              "leachatelow"                    
     [5] "leachatemid"                     "leachatehigh"                   
     [7] "stageprawnchip.leachatelow"      "stageearlygastrula.leachatelow" 
     [9] "stageprawnchip.leachatemid"      "stageearlygastrula.leachatemid" 
    [11] "stageprawnchip.leachatehigh"     "stageearlygastrula.leachatehigh"

``` r
# cleavage 
res_cleavage_low.vs.control <- results(dds_wald, name = "leachatelow", alpha = 0.05)
res_cleavage_mid.vs.control  <- results(dds_wald, name = "leachatemid", alpha = 0.05)
res_cleavage_high.vs.control <- results(dds_wald, name = "leachatehigh", alpha = 0.05)

# prawnchip
res_prawnchip_low.vs.control <- results(dds_wald, contrast = list(c("leachatelow", "stageprawnchip.leachatelow")), alpha = 0.05)
res_prawnchip_mid.vs.control  <- results(dds_wald, contrast = list(c("leachatemid",  "stageprawnchip.leachatemid")), alpha = 0.05)
res_prawnchip_high.vs.control <- results(dds_wald, contrast = list(c("leachatehigh", "stageprawnchip.leachatehigh")), alpha = 0.05)

# earlygastrula
res_earlygastrula_low.vs.control <- results(dds_wald, contrast = list(c("leachatelow", "stageearlygastrula.leachatelow")), alpha = 0.05)
res_earlygastrula_mid.vs.control  <- results(dds_wald, contrast = list(c("leachatemid",  "stageearlygastrula.leachatemid")), alpha = 0.05)
res_earlygastrula_high.vs.control <- results(dds_wald, contrast = list(c("leachatehigh", "stageearlygastrula.leachatehigh")), alpha = 0.05)
```

### Cleavage

#### high

``` r
summary(res_cleavage_high.vs.control)
```


    out of 12565 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0 (up)       : 5, 0.04%
    LFC < 0 (down)     : 13, 0.1%
    outliers [1]       : 0, 0%
    low counts [2]     : 0, 0%
    (mean count < 23)
    [1] see 'cooksCutoff' argument of ?results
    [2] see 'independentFiltering' argument of ?results

> [!NOTE]
>
> 18 genes were found to be significantly differentially expressed
> (adjusted p value \< 0.05) at the cleavage stage in the high leachate
> treatment compared to the control. 13 LFC \< 0. These 13 transcripts
> were all found in less abundance compared to the controls.

``` r
sig_cleavage_high.vs.control <- res_cleavage_high.vs.control %>%
  as.data.frame() %>%
  filter(padj < 0.05) %>% 
  rownames_to_column(var = "gene_id")

write_csv(sig_cleavage_high.vs.control, file.path(output, "sig_cleavage_high.vs.control.csv"))
```

#### mid

mid leachate vs. control (mvc) for cleavage

``` r
summary(res_cleavage_mid.vs.control)
```


    out of 12565 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0 (up)       : 77, 0.61%
    LFC < 0 (down)     : 0, 0%
    outliers [1]       : 0, 0%
    low counts [2]     : 0, 0%
    (mean count < 23)
    [1] see 'cooksCutoff' argument of ?results
    [2] see 'independentFiltering' argument of ?results

> [!NOTE]
>
> 77 genes were found to be significantly differentially expressed
> (adjusted p value \< 0.05) at the cleavage stage in the mid leachate
> treatment compared to the control. 77 LFC \> 0. These 77 genes were
> all found in more abundance compared to the controls.

``` r
sig_cleavage_mid.vs.control <- res_cleavage_mid.vs.control %>%
  as.data.frame() %>%
  filter(padj < 0.05) %>% 
  rownames_to_column(var = "gene_id")

write_csv(sig_cleavage_mid.vs.control, file.path(output, "sig_cleavage_mid.vs.control.csv"))
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


    out of 12565 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0 (up)       : 14, 0.11%
    LFC < 0 (down)     : 25, 0.2%
    outliers [1]       : 0, 0%
    low counts [2]     : 0, 0%
    (mean count < 23)
    [1] see 'cooksCutoff' argument of ?results
    [2] see 'independentFiltering' argument of ?results

> [!NOTE]
>
> 39 genes were found to be significantly differentially expressed
> (adjusted p value \< 0.05) at the cleavage stage in the low leachate
> treatment compared to the control. 14 LFC \> 0 & 25 LFC \< 0.

``` r
sig_cleavage_low.vs.control <- res_cleavage_low.vs.control %>%
  as.data.frame() %>%
  filter(padj < 0.05) %>% 
  rownames_to_column(var = "gene_id")

write_csv(sig_cleavage_low.vs.control, file.path(output, "sig_cleavage_low.vs.control.csv"))
```

### Prawnchip

#### high

for prawnchip

``` r
summary(res_prawnchip_high.vs.control)
```


    out of 12565 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0 (up)       : 0, 0%
    LFC < 0 (down)     : 0, 0%
    outliers [1]       : 0, 0%
    low counts [2]     : 0, 0%
    (mean count < 23)
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


    out of 12565 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0 (up)       : 0, 0%
    LFC < 0 (down)     : 0, 0%
    outliers [1]       : 0, 0%
    low counts [2]     : 0, 0%
    (mean count < 23)
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


    out of 12565 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0 (up)       : 0, 0%
    LFC < 0 (down)     : 0, 0%
    outliers [1]       : 0, 0%
    low counts [2]     : 0, 0%
    (mean count < 23)
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


    out of 12565 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0 (up)       : 22, 0.18%
    LFC < 0 (down)     : 0, 0%
    outliers [1]       : 0, 0%
    low counts [2]     : 10231, 81%
    (mean count < 1572)
    [1] see 'cooksCutoff' argument of ?results
    [2] see 'independentFiltering' argument of ?results

> [!NOTE]
>
> 22 genes differentially expressed (adjusted p value \< 0.05) at the
> early gastrula stage in the high leachate treatment when compared to
> the control. All 22 of these genes had more abundant transcripts in
> treated samples compared to controls.

``` r
sig_earlygastrula_high.vs.control <- res_earlygastrula_high.vs.control %>%
  as.data.frame() %>%
  filter(padj < 0.05) %>% 
  rownames_to_column(var = "gene_id")

write_csv(sig_earlygastrula_high.vs.control, file.path(output, "sig_earlygastrula_high.vs.control.csv"))
```

#### mid

``` r
summary(res_earlygastrula_mid.vs.control)
```


    out of 12565 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0 (up)       : 0, 0%
    LFC < 0 (down)     : 0, 0%
    outliers [1]       : 0, 0%
    low counts [2]     : 0, 0%
    (mean count < 23)
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


    out of 12565 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0 (up)       : 0, 0%
    LFC < 0 (down)     : 0, 0%
    outliers [1]       : 0, 0%
    low counts [2]     : 0, 0%
    (mean count < 23)
    [1] see 'cooksCutoff' argument of ?results
    [2] see 'independentFiltering' argument of ?results

> [!NOTE]
>
> 0 genes differentially expressed (adjusted p value \< 0.05) at the
> early gastrula stage in the low leachate treatment when compared to
> the control.

# Summary

So we have 11 different ‘significant gene’ lists. Each list indicates
that leachate had a significant effect somewhere:

LRT - `sig_cat_LRT`, 44 DEGs, this list are genes that had an overall
impact by leachate interacting with stage (full vs reduced, testing
interaction term as a whole)

Interaction effects - `sig_prawnchip_high.int`, 1 DEG -
`sig_prawnchip_mid.int`, 5 DEGs - `sig_prawnchip_low.int`, 12 DEGs -
`sig_earlygastrula_high.int`, 3 DEGs - `sig_earlygastrula_mid.int`, 11
DEGs - `sig_earlygastrula_low.int`, 11 DEGs

Stage x leachate effects - `sig_cleavage_low.vs.control`, 39 DEGs -
`sig_cleavage_mid.vs.control`, 77 DEGs - `sig_cleavage_high.vs.control`,
18 DEGs - `sig_earlygastrula_high.vs.control`, 22 DEGs

``` r
# 1) Put all of your sig data frames into a *named* list
sig_lists <- list(
  LRT                         = sig_cat_LRT,
  prawnchip_high.int          = sig_prawnchip_high.int,
  prawnchip_mid.int           = sig_prawnchip_mid.int,
  prawnchip_low.int           = sig_prawnchip_low.int,
  earlygastrula_high.int      = sig_earlygastrula_high.int,
  earlygastrula_mid.int       = sig_earlygastrula_mid.int,
  earlygastrula_low.int       = sig_earlygastrula_low.int,
  cleavage_low.vs.control     = sig_cleavage_low.vs.control,
  cleavage_mid.vs.control     = sig_cleavage_mid.vs.control,
  cleavage_high.vs.control    = sig_cleavage_high.vs.control,
  earlygastrula_high.vs.control = sig_earlygastrula_high.vs.control
)

# 2) Smoosh: add a 'list' column indicating the source and stack them
sig_summary <- bind_rows(sig_lists, .id = "list") %>%
  relocate(list, gene_id)

write_csv(sig_summary, file.path(output, "sig_summary.csv"))

# Quick sanity check of counts you reported
sig_counts <- sig_summary %>% count(list, name = "n_genes")
sig_counts
```

                                list n_genes
    1                            LRT      44
    2       cleavage_high.vs.control      18
    3        cleavage_low.vs.control      39
    4        cleavage_mid.vs.control      77
    5         earlygastrula_high.int       3
    6  earlygastrula_high.vs.control      22
    7          earlygastrula_low.int      11
    8          earlygastrula_mid.int      11
    9             prawnchip_high.int       1
    10             prawnchip_low.int      12
    11             prawnchip_mid.int       5

# Next steps

- I’d simply like an easy absence/presence heatmap that shows all the
  genes, and which lists they are significant in.
- I’d then like to show on the heatmap the log2FC, and add astericks for
  significance levels p\<0.05, p\<0.01, p\<0.001.
- I’d like gene count plots visually summarizing each list type (LRT,
  Interaction, Stagexleachate effects) in tiny tidy multiples
- 
