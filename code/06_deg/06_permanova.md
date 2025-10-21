# PERMANOVA of RNA-seq gene counts
Sarah Tanja
2024-12-02

- [<span class="toc-section-number">1</span> Background](#background)
  - [<span class="toc-section-number">1.1</span> Inputs](#inputs)
  - [<span class="toc-section-number">1.2</span> Outputs](#outputs)
  - [<span class="toc-section-number">1.3</span> Load
    packages](#load-packages)
  - [<span class="toc-section-number">1.4</span> Load
    inputs](#load-inputs)
- [<span class="toc-section-number">2</span> Make DESeqData
  Object](#make-deseqdata-object)
- [<span class="toc-section-number">3</span> Full
  PERMANOVA](#full-permanova)
  - [<span class="toc-section-number">3.1</span> Estimate size
    factors](#estimate-size-factors)
  - [<span class="toc-section-number">3.2</span> Apply a variance
    stabilizing transformation
    (VST)](#apply-a-variance-stabilizing-transformation-vst)
  - [<span class="toc-section-number">3.3</span> Compute distance
    between samples](#compute-distance-between-samples)
  - [<span class="toc-section-number">3.4</span> Run `adonis2` to build
    PERMANOVA model](#run-adonis2-to-build-permanova-model)
- [<span class="toc-section-number">4</span> DEG
  PERMANOVA](#deg-permanova)
  - [<span class="toc-section-number">4.1</span> Subset DESeqDataSet by
    significant gene list
    `sig_all`](#subset-deseqdataset-by-significant-gene-list-sig_all)
  - [<span class="toc-section-number">4.2</span> Estimate size
    factors](#estimate-size-factors-1)
  - [<span class="toc-section-number">4.3</span> Apply variance
    stabilizing
    transformation](#apply-variance-stabilizing-transformation)
  - [<span class="toc-section-number">4.4</span> Compute distances
    between samples](#compute-distances-between-samples)
  - [<span class="toc-section-number">4.5</span> Run adonis2 to build
    PERMANOVA model](#run-adonis2-to-build-permanova-model-1)

# Background

## Inputs

- metadata `metadata_or`
- filtered gene count matrix `gcm_filtor`
- list of all significantly expressed genes `sig_all`

## Outputs

## Load packages

``` r
library(DESeq2)
library(vegan)
library(tidyverse)
```

## Load inputs

``` r
sig_all <- read_csv("../../output/06_deg/results/sig_all.csv")
```

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

``` r
# Load metadata (outliers removed)
metadata_or <- read_csv("../../metadata/metadata_or.csv")

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

# Make DESeqData Object

``` r
dds_cat <- DESeqDataSetFromMatrix(countData = gcm_filtor,
                              colData = metadata_or,
                              design = formula( ~0 + stage + leachate + stage:leachate))
```

# Full PERMANOVA

## Estimate size factors

check that they are all less than 4 (should return TRUE)

``` r
dds_cat <- estimateSizeFactors(dds_cat)
hist(sizeFactors(dds_cat))
```

<img src="06_permanova_files/figure-commonmark/unnamed-chunk-6-1.jpeg"
width="4200" />

``` r
all(sizeFactors(dds_cat))<4
```

    [1] TRUE

## Apply a variance stabilizing transformation (VST)

We apply a variance stabilizing transformation to minimize effects of
small counts and normalize for library size

Note: t() transposes. With DESeq2’s VST:

t(assay(vst(dds_cat, blind = TRUE))) flips it to samples × genes.

Why it matters:

Functions like vegan::vegdist() and adonis2() expect rows = samples,
columns = features.

Use the transposed version for PERMANOVA: rows (samples) will match your
meta data frame.

``` r
# extract variance stabilized count matrix
vst_mat <- assay(vst(dds_cat, blind = TRUE))

# transpose
vst_t_mat <- t(vst_mat)
str(vst_t_mat)
```

     num [1:61, 1:12565] 7.14 9.66 8.23 7.42 8.33 ...
     - attr(*, "dimnames")=List of 2
      ..$ : chr [1:61] "101112C14" "101112C4" "101112C9" "101112H14" ...
      ..$ : chr [1:12565] "Montipora_capitata_HIv3___RNAseq.g4581.t1" "Montipora_capitata_HIv3___RNAseq.g4582.t1" "Montipora_capitata_HIv3___RNAseq.g4583.t1" "Montipora_capitata_HIv3___RNAseq.g4584.t1" ...

## Compute distance between samples

``` r
eu_dist <- vegdist((vst_t_mat), method = "eu")
```

## Run `adonis2` to build PERMANOVA model

> [!NOTE]
>
> Because adonis2() uses random permutations to get p-values. Without
> fixing the seed, each run draws a different set/order of permutations,
> which results in tiny differences in F-stats/p-values. By adding
> `set.seed(1)` the permutation stream is reproducible so you (and
> future you) get the same numbers.

``` r
adonis2(eu_dist ~ leachate * stage, 
        data = metadata_or,     
        permutations = 999,
        by = "terms")
```

    Permutation test for adonis under reduced model
    Terms added sequentially (first to last)
    Permutation: free
    Number of permutations: 999

    adonis2(formula = eu_dist ~ leachate * stage, data = metadata_or, permutations = 999, by = "terms")
                   Df SumOfSqs      R2       F Pr(>F)    
    leachate        3     9640 0.02234  1.0947  0.337    
    stage           2   264333 0.61263 45.0264  0.001 ***
    leachate:stage  6    13666 0.03167  0.7760  0.715    
    Residual       49   143831 0.33335                   
    Total          60   431470 1.00000                   
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
adonis2(eu_dist ~ leachate * stage, 
        data = metadata_or,
        permutations = 999,
        strata = metadata_or$parents, # <-- control for parental cross
        by = "terms") 
```

    Permutation test for adonis under reduced model
    Terms added sequentially (first to last)
    Blocks:  strata 
    Permutation: free
    Number of permutations: 999

    adonis2(formula = eu_dist ~ leachate * stage, data = metadata_or, permutations = 999, by = "terms", strata = metadata_or$parents)
                   Df SumOfSqs      R2       F Pr(>F)    
    leachate        3     9640 0.02234  1.0947  0.383    
    stage           2   264333 0.61263 45.0264  0.001 ***
    leachate:stage  6    13666 0.03167  0.7760  0.701    
    Residual       49   143831 0.33335                   
    Total          60   431470 1.00000                   
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

> [!NOTE]
>
> A PERMANOVA with 999 permutations (permutations restricted within
> parental cross) revealed that developmental stage explained the
> majority of variation in multivariate gene-expression profiles
> ($F_{2,49} = 45.0$, $R^2 = 0.613$, $p = 0.001$). Neither leachate
> concentration ($F_{3,49} = 1.09$, $R^2 = 0.022$, $p = 0.375$) nor the
> leachate,×,stage interaction ($F_{6,49} = 0.78$, $R^2 = 0.032$,
> $p = 0.697$) had a significant effect. Overall, the model explained
> 66.7% of the variance in gene-expression dissimilarity, with stage
> accounting for nearly all of the explained variation.
>
> These results indicate that transcriptomic profiles differed strongly
> among developmental stages but not across leachate concentrations or
> their interaction, suggesting that developmental processes dominate
> gene-expression variation in this dataset.

Why does the permanova say `leachate` and `leachate:stage` is
insignificant but `DESeq2` results return significantly differentially
expressed genes for those terms?

**PERMANOVA** is a **single** multivariate test on a distance among
**samples**

It asks whether group **centroids** differ (and a bit of spread,
depending on dispersion). If `leachate` shifts are modest,
bidirectional, or limited to a subset of genes, the overall centroid
shift can be tiny → non-sig.

**Dominant developmental stage signal**

Your term table shows **`stage`** explains ~61% of variance; that can
swamp smaller `leachate` signals. After blocking by parents, the
residual degrees of freedom and variability left to detect `leachate`
are limited.

Your DESeq2 results might be **stage-specific** contrasts (e.g.,
leachate within stage), while your PERMANOVA tested a **single global**
leachate and interaction term. The global test can be ns even if **some
stages** have effects.

# DEG PERMANOVA

``` r
head(rownames(dds_cat))
```

    [1] "Montipora_capitata_HIv3___RNAseq.g4581.t1"
    [2] "Montipora_capitata_HIv3___RNAseq.g4582.t1"
    [3] "Montipora_capitata_HIv3___RNAseq.g4583.t1"
    [4] "Montipora_capitata_HIv3___RNAseq.g4584.t1"
    [5] "Montipora_capitata_HIv3___RNAseq.g4585.t1"
    [6] "Montipora_capitata_HIv3___RNAseq.g4586.t1"

``` r
 head(sig_all)
```

    # A tibble: 6 × 1
      gene_id                                  
      <chr>                                    
    1 Montipora_capitata_HIv3___TS.g26493.t2a  
    2 Montipora_capitata_HIv3___RNAseq.g5198.t1
    3 Montipora_capitata_HIv3___TS.g28558.t3a  
    4 Montipora_capitata_HIv3___RNAseq.g6929.t1
    5 Montipora_capitata_HIv3___RNAseq.g7017.t1
    6 Montipora_capitata_HIv3___TS.g42133.t1   

## Subset DESeqDataSet by significant gene list `sig_all`

``` r
dds_mat <- counts(dds_cat) 

keep <-  rownames(dds_mat) %in% sig_all$gene_id

dds_deg <- dds_cat[keep, , drop = FALSE]  
```

## Estimate size factors

``` r
dds_deg <- estimateSizeFactors(dds_deg)
hist(sizeFactors(dds_deg))
```

<img src="06_permanova_files/figure-commonmark/unnamed-chunk-14-1.jpeg"
width="4200" />

``` r
all(sizeFactors(dds_cat))<4
```

    [1] TRUE

## Apply variance stabilizing transformation

``` r
# extract variance stabilized count matrix
vst_deg <- assay(varianceStabilizingTransformation(dds_deg, blind = TRUE, fitType = "local"))

# transpose
vst_t_deg <- t(vst_deg)
str(vst_t_deg)
```

     num [1:61, 1:151] 8.05 -6.43 4.32 7.69 1.49 ...
     - attr(*, "dimnames")=List of 2
      ..$ : chr [1:61] "101112C14" "101112C4" "101112C9" "101112H14" ...
      ..$ : chr [1:151] "Montipora_capitata_HIv3___TS.g26493.t2a" "Montipora_capitata_HIv3___RNAseq.g5198.t1" "Montipora_capitata_HIv3___RNAseq.g5394.t1" "Montipora_capitata_HIv3___TS.g28558.t3a" ...

## Compute distances between samples

``` r
eu_dist_deg <- vegdist((vst_t_deg), method = "eu")
```

## Run adonis2 to build PERMANOVA model

``` r
adonis2(eu_dist_deg ~ leachate * stage, 
        data = metadata_or,
        permutations = 999,
        strata = metadata_or$parents, # <-- control for parental cross
        by = "terms") 
```

    Permutation test for adonis under reduced model
    Terms added sequentially (first to last)
    Blocks:  strata 
    Permutation: free
    Number of permutations: 999

    adonis2(formula = eu_dist_deg ~ leachate * stage, data = metadata_or, permutations = 999, by = "terms", strata = metadata_or$parents)
                   Df SumOfSqs      R2       F Pr(>F)    
    leachate        3     4453 0.03930  2.9581  0.035 *  
    stage           2    77096 0.68046 76.8251  0.001 ***
    leachate:stage  6     7164 0.06323  2.3796  0.035 *  
    Residual       49    24586 0.21700                   
    Total          60   113299 1.00000                   
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

> [!NOTE]
>
> Using PERMANOVA on Euclidean distances calculated from the 151
> differentially expressed genes (DEGs) (n = 61 samples; 999
> permutations; stratified by parental cross), we detected a strong
> effect of developmental stage and smaller but significant effects of
> leachate and their interaction. Stage explained R² = 0.680 of the
> variance and was highly significant (Df = 2, F = 76.83, p = 0.001).
> Leachate explained R² = 0.039 (Df = 3, F = 2.96, p = 0.034), and the
> leachate × stage interaction explained R² = 0.063 (Df = 6, F = 2.38, p
> = 0.023). Residual variation accounted for R² = 0.217.
