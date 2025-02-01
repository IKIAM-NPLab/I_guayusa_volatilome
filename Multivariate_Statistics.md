Effect of age and light on the volatile profile of *Ilex guayusa*
leaves - Statistics analysis
================
Jefferson Pastuña
2024-04-17

- <a href="#introduction" id="toc-introduction">Introduction</a>
- <a href="#data-preprocessing" id="toc-data-preprocessing">Data
  preprocessing</a>
  - <a href="#before-to-start" id="toc-before-to-start">Before to start</a>
  - <a href="#notame-workflow" id="toc-notame-workflow">Notame workflow</a>
  - <a href="#preprocessing" id="toc-preprocessing">Preprocessing</a>
- <a href="#univariate-statistic" id="toc-univariate-statistic">Univariate
  statistic</a>
  - <a href="#volcano-plot" id="toc-volcano-plot">Volcano plot</a>
    - <a href="#subqc-volcano-plot" id="toc-subqc-volcano-plot">SubQC volcano
      plot</a>
    - <a href="#sample-volcano-plot" id="toc-sample-volcano-plot">Sample
      volcano plot</a>
  - <a href="#tukey-test" id="toc-tukey-test">Tukey test</a>
    - <a href="#subqc-tukey-test" id="toc-subqc-tukey-test">SubQC Tukey
      test</a>
    - <a href="#sample-tukey-test" id="toc-sample-tukey-test">Sample Tukey
      test</a>
- <a href="#multivariate-statistic"
  id="toc-multivariate-statistic">Multivariate statistic</a>
  - <a href="#principal-component-analysis-pca"
    id="toc-principal-component-analysis-pca">Principal Component Analysis
    (PCA)</a>
    - <a href="#light-factor-pca" id="toc-light-factor-pca">Light factor
      PCA</a>
    - <a href="#age-factor-pca" id="toc-age-factor-pca">Age factor PCA</a>
    - <a href="#location-factor-pca" id="toc-location-factor-pca">Location
      factor PCA</a>
    - <a href="#alto-pano-location-pca" id="toc-alto-pano-location-pca">Alto
      Pano location PCA</a>
    - <a href="#alto-tena-location-pca" id="toc-alto-tena-location-pca">Alto
      Tena location PCA</a>
    - <a href="#talag-location-pca" id="toc-talag-location-pca">Talag location
      PCA</a>
  - <a href="#heatmap-and-hca" id="toc-heatmap-and-hca">Heatmap and HCA</a>
    - <a href="#heatmap-and-hca-of-age"
      id="toc-heatmap-and-hca-of-age">Heatmap and HCA of age</a>
    - <a href="#heatmap-and-hca-of-sample"
      id="toc-heatmap-and-hca-of-sample">Heatmap and HCA of sample</a>

# Introduction

This R Script aims to record the procedure for analyzing the volatile
profile of *Ilex guayusa* leaves under different age and light
conditions. Each step has a brief explanation, as well as code and
graphics.

The data preprocessing workflow used was taken from [“notame”: Workflow
for Non-Targeted LC–MS Metabolic
Profiling](https://doi.org/10.3390/metabo10040135), Which offers a wide
variety of functions for perform metabolomic profile analysis.

# Data preprocessing

## Before to start

The “notame” package accepts as input a feature table that can be
obtained through software such as MZmine, MS-DIAL, among others. In this
case, the feature table was obtained with the help of MZmine. The
(\*.csv) file exported from MZmine was fixed to get the final feature
table input according to the “notame” package format.

Modifications to the raw (\*.csv) file can be summarized by adding and
renaming columns. The added columns “Column” and “Ion Mode” allow for
the analysis of samples with different types of columns and different
ionization modes, respectively. Also, the cells corresponding to mass
and retention time must be renamed so the “notame” package can detect
and process them.

## Notame workflow

The “notame” package and other dependency packages were installed as a
first step for the analysis.

``` r
# Notame package installation
#if (!requireNamespace("devtools", quietly = TRUE)) {
#  install.packages("devtools")
#}
#devtools::install_github("antonvsdata/notame", ref = "v0.3.1")

# Notame library call
library(notame)

# Dependency packages installation
install_dependencies
```

    ## function (preprocessing = TRUE, extra = FALSE, batch_corr = FALSE, 
    ##     misc = FALSE, ...) 
    ## {
    ##     core_cran <- c("BiocManager", "cowplot", "missForest", "openxlsx", 
    ##         "randomForest", "RColorBrewer", "Rtsne")
    ##     core_bioconductor <- "pcaMethods"
    ##     extra_cran <- c("car", "doParallel", "devEMF", "ggbeeswarm", 
    ##         "ggdendro", "ggrepel", "ggtext", "Hmisc", "hexbin", "igraph", 
    ##         "lme4", "lmerTest", "MuMIn", "PERMANOVA", "PK", "rmcorr")
    ##     extra_bioconductor <- c("mixOmics", "supraHex")
    ##     extra_gitlab <- "CarlBrunius/MUVR"
    ##     batch_cran <- "fpc"
    ##     batch_bioconductor <- "RUVSeq"
    ##     batch_github <- NULL
    ##     batch_gitlab <- "CarlBrunius/batchCorr"
    ##     misc_cran <- c("knitr", "rmarkdown", "testthat")
    ##     if (preprocessing) {
    ##         install_helper(cran = core_cran, bioconductor = core_bioconductor, 
    ##             ...)
    ##     }
    ##     if (extra) {
    ##         install_helper(cran = extra_cran, bioconductor = extra_bioconductor, 
    ##             gitlab = extra_gitlab, ...)
    ##     }
    ##     if (batch_corr) {
    ##         install_helper(cran = batch_cran, bioconductor = batch_bioconductor, 
    ##             github = batch_github, gitlab = batch_gitlab, ...)
    ##     }
    ##     if (misc) {
    ##         install_helper(cran = misc_cran, ...)
    ##     }
    ## }
    ## <bytecode: 0x000000000fc40de0>
    ## <environment: namespace:notame>

Then, a primary path and a log system were added to have a record of
each process executed.

``` r
# Main path
ppath <- "../I_guayusa_volatilome/"
# Log system
init_log(log_file = paste0(ppath, "Result/notame_Result/GCMS_log.txt"))
```

    ## INFO [2025-01-31 19:24:22] Starting logging

Next, the MZmine feature list in “notame” format was loaded.

``` r
data <- read_from_excel(file = "Data/Data_to_notame/MZmine_Feature_List_2Notame.xlsx", sheet = 3, 
                        corner_row = 8, corner_column = "N", 
                        split_by = c("Column", "Ion mode"))
```

    ## INFO [2025-01-31 19:24:24] Corner detected correctly at row 8, column N
    ## INFO [2025-01-31 19:24:24] 
    ## Extracting sample information from rows 1 to 8 and columns O to DC
    ## INFO [2025-01-31 19:24:24] Replacing spaces in sample information column names with underscores (_)
    ## INFO [2025-01-31 19:24:24] Naming the last column of sample information "Datafile"
    ## INFO [2025-01-31 19:24:24] 
    ## Extracting feature information from rows 9 to 5104 and columns A to N
    ## INFO [2025-01-31 19:24:24] Creating Split column from Column, Ion mode
    ## INFO [2025-01-31 19:24:24] Feature_ID column not found, creating feature IDs
    ## INFO [2025-01-31 19:24:24] Identified m/z column mass and retention time column rt
    ## INFO [2025-01-31 19:24:24] Identified m/z column mass and retention time column rt
    ## INFO [2025-01-31 19:24:24] Creating feature IDs from Split, m/z and retention time
    ## INFO [2025-01-31 19:24:24] Replacing dots (.) in feature information column names with underscores (_)
    ## INFO [2025-01-31 19:24:24] 
    ## Extracting feature abundances from rows 9 to 5104 and columns O to DC
    ## INFO [2025-01-31 19:24:24] 
    ## Checking sample information
    ## INFO [2025-01-31 19:24:24] QC column generated from rows containing 'QC'
    ## INFO [2025-01-31 19:24:24] Sample ID autogenerated from injection orders and prefix ID_
    ## INFO [2025-01-31 19:24:24] Checking that feature abundances only contain numeric values
    ## INFO [2025-01-31 19:24:25] 
    ## Checking feature information
    ## INFO [2025-01-31 19:24:25] Checking that feature IDs are unique and not stored as numbers
    ## INFO [2025-01-31 19:24:25] Checking that m/z and retention time values are reasonable
    ## INFO [2025-01-31 19:24:25] Identified m/z column mass and retention time column rt
    ## INFO [2025-01-31 19:24:25] Identified m/z column mass and retention time column rt

Once the data was loaded, the next step was to create a MetaboSet to
work with R objects from now on.

``` r
modes <- construct_metabosets(exprs = data$exprs, 
                              pheno_data = data$pheno_data, 
                              feature_data = data$feature_data,
                              group_col = "Group")
```

    ## Initializing the object(s) with unflagged features
    ## INFO [2025-01-31 19:24:25] 
    ## Checking feature information
    ## INFO [2025-01-31 19:24:25] Checking that feature IDs are unique and not stored as numbers
    ## INFO [2025-01-31 19:24:25] Checking that feature abundances only contain numeric values
    ## INFO [2025-01-31 19:24:25] Setting row and column names of exprs based on feature and pheno data

## Preprocessing

The first step of the preprocessing is to change the features with a
value equal to 0 to NA.

``` r
# Data extraction
mode <- modes$Rtx5MS_EI
# Change 0 value to NA
mode <- mark_nas(mode, value = 0)
```

Then, features with low detection rates are flagged and can be ignored
or removed in subsequent analysis. The “notame” package uses two
criteria to flag these features: the feature presence in a percentage of
QC injections and the feature presence in a percentage within a sample
group or class.

``` r
# Low detection rate
mode <- flag_detection(mode, qc_limit = 10/14, group_limit = 2/3)
```

    ## INFO [2025-01-31 19:24:27] 
    ## 0% of features flagged for low detection rate

The following preprocessing step is drift correction, which is applied
using smoothed cubic spline regression. In addition, features with low
quality (features out of the “RSD_r \< 0.3 & D_ratio_r \< 0.6”
condition) are added to the “Flag” column. The value of RSD \< 0.3 for
considering a low-quality feature was taken from [Procedures for
large-scale metabolic profiling of serum and plasma using gas
chromatography and liquid chromatography coupled to mass
spectrometry](https://doi.org/10.1038/nprot.2011.335).

``` r
# Drift correction
corrected <- correct_drift(mode)
```

    ## INFO [2025-01-31 19:24:27] 
    ## Starting drift correction at 2025-01-31 19:24:27.960347
    ## INFO [2025-01-31 19:25:22] Drift correction performed at 2025-01-31 19:25:22.604473
    ## INFO [2025-01-31 19:25:29] Inspecting drift correction results 2025-01-31 19:25:29.562871
    ## INFO [2025-01-31 19:26:01] Drift correction results inspected at 2025-01-31 19:26:01.593703
    ## INFO [2025-01-31 19:26:01] 
    ## Drift correction results inspected, report:
    ## Drift_corrected: 100%,  Missing_QCS: 0%

``` r
# Flag low quality features
corrected <- flag_quality(corrected,
                          condition = "RSD_r < 0.3 & D_ratio_r < 0.6")
```

    ## INFO [2025-01-31 19:26:01] 
    ## 12% of features flagged for low quality

Then, we can visualize the data after drift correction.

``` r
# Boxplot
corr_bp <- plot_sample_boxplots(corrected,
                                order_by = "QC",
                                fill_by = "QC")
# PCA
corr_pca <- plot_pca(corrected,
                     center = TRUE,
                     shape = "QC",
                     color = "QC")
# Package to plots visualization in a same windows
#if (!requireNamespace("devtools", quietly = TRUE)) {
#  install.packages("devtools")
#}
#devtools::install_github("thomasp85/patchwork")
library(patchwork)
# Plot
corr_pca + corr_bp
```

![](Multivariate_Statistics_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

Contaminant peaks based on the process blank sample will be removed.

``` r
# Removal of contaminants
corrected_no_blank <- flag_contaminants(corrected,
                                        blank_col = "Group",
                                        blank_label = "Blank",
                                        flag_thresh = 0.39,
                                        flag_label = "Contaminant")
```

    ## INFO [2025-01-31 19:26:06] 
    ## 42% of features flagged as contaminants

``` r
# Removal blank group from dataset
corrected_no_blank <- corrected_no_blank[, corrected_no_blank$Group != "Blank"]
pData(corrected_no_blank) <- droplevels(pData(corrected_no_blank))
# Exporting data to clean identified features
#write_to_excel(corrected_no_blank,
#               "Result/notame_Result/MZmine_corrected_no_blank.xlsx")
```

We can inspect data after removing the contaminant features.

``` r
# Boxplot
no_blank_bp <- plot_sample_boxplots(corrected_no_blank,
                                    order_by = "QC",
                                    fill_by = "QC")
# PCA
no_blank_pca <- plot_pca(corrected_no_blank,
                         center = TRUE,
                         shape = "QC",
                         color = "QC")
# Plot
no_blank_pca + no_blank_bp
```

![](Multivariate_Statistics_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

The next step is feature clustering. This step helps us reduce the
number of features of the same molecule that were split due to a 70 eV
electron ionization (EI) environment.

``` r
# Clustering correlated features
clustered <- cluster_features(corrected_no_blank,
                              rt_window = 1/180,
                              corr_thresh = 0.95,
                              d_thresh = 0.80,
                              #plotting = TRUE,
                              #prefix = paste0(ppath, "Result/notame_Result/Cluster/")
                              )
```

    ## INFO [2025-01-31 19:26:09] Identified m/z column mass and retention time column rt
    ## INFO [2025-01-31 19:26:09] 
    ## Starting feature clustering at 2025-01-31 19:26:09.94118
    ## INFO [2025-01-31 19:26:09] Finding connections between features in Rtx5MS_EI
    ## [1] 100
    ## [1] 200
    ## [1] 300
    ## [1] 400
    ## [1] 500
    ## [1] 600
    ## [1] 700
    ## [1] 800
    ## [1] 900
    ## [1] 1000
    ## [1] 1100
    ## [1] 1200
    ## [1] 1300
    ## [1] 1400
    ## [1] 1500
    ## [1] 1600
    ## [1] 1700
    ## [1] 1800
    ## [1] 1900
    ## [1] 2000
    ## [1] 2100
    ## [1] 2200
    ## [1] 2300
    ## INFO [2025-01-31 19:29:55] Found 654 connections in Rtx5MS_EI
    ## INFO [2025-01-31 19:29:55] Found 654 connections
    ## 245 components found
    ## 
    ## Component 100 / 245 
    ## Component 200 / 245 
    ## 101 components found
    ## 
    ## Component 100 / 101 
    ## 4 components found
    ## 
    ## INFO [2025-01-31 19:29:56] Found 267 clusters of 2 or more features, clustering finished at 2025-01-31 19:29:56.104116

``` r
# Extracting representative feature of each cluster
compressed <- compress_clusters(clustered)
```

    ## INFO [2025-01-31 19:29:56] Clusters compressed, left with 1939 features

``` r
# Exporting data to inspect identified features
#write_to_excel(compressed,
#               "Result/notame_Result/MZmine_compressed.xlsx")
```

We can inspect the data using the PCA plot after the clustering
algorithm execution.

``` r
# Boxplot
clust_bp <- plot_sample_boxplots(compressed,
                                 order_by = "QC",
                                 fill_by = "QC")
# PCA
clust_pca <- plot_pca(compressed,
                      center = TRUE,
                      shape = "QC",
                      color = "QC")
# Plot
clust_pca + clust_bp
```

![](Multivariate_Statistics_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

For downstream statistical analysis, we used probabilistic quotient
normalization (PQN) with random forest missing value imputation and
generalized logarithm (glog) transformation method, which were taken
from [Non-targeted UHPLC-MS metabolomic data processing methods: a
comparative investigation of normalisation, missing value imputation,
transformation and scaling](https://doi.org/10.1007/s11306-016-1030-9).

The code below imputes the data using a random forest algorithm.

``` r
# Impute missing values using random forest
# To clean data
set.seed(1010)
imputed <- impute_rf(compressed)
```

    ## INFO [2025-01-31 19:29:59] 
    ## Starting random forest imputation at 2025-01-31 19:29:59.175292
    ## INFO [2025-01-31 19:30:52] Out-of-bag error in random forest imputation: 0.324
    ## INFO [2025-01-31 19:30:52] Random forest imputation finished at 2025-01-31 19:30:52.921366

``` r
# To all data
imputed <- impute_rf(imputed, all_features = TRUE)
```

    ## INFO [2025-01-31 19:30:52] 
    ## Starting random forest imputation at 2025-01-31 19:30:52.925366
    ## INFO [2025-01-31 19:30:53] Out-of-bag error in random forest imputation: 0
    ## INFO [2025-01-31 19:30:53] Random forest imputation finished at 2025-01-31 19:30:53.009371

We can inspect the data with the PCA plot after data imputation.

``` r
# Boxplot
imp_bp <- plot_sample_boxplots(imputed,
                               order_by = "QC",
                               fill_by = "QC")
# PCA
imp_pca <- plot_pca(imputed,
                    center = TRUE,
                    shape = "QC",
                    color = "QC")
# Plot
imp_pca + imp_bp
```

![](Multivariate_Statistics_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

After data imputation, the data was normalized by probabilistic quotient
normalization (PQN).

``` r
# Probabilistic quotient normalization
pqn_set <- pqn_normalization(imputed,
                             ref = c("qc", "all"),
                             method = c("median", "mean"),
                             all_features = FALSE)
```

    ## INFO [2025-01-31 19:30:55] Starting PQN normalization
    ## INFO [2025-01-31 19:30:55] Using median of qc samples as reference spectrum

We can inspect the data with the PCA plot after data normalization.

``` r
# Boxplot
pqn_bp <- plot_sample_boxplots(pqn_set,
                               order_by = "QC",
                               fill_by = "QC")
# PCA
pqn_pca <- plot_pca(pqn_set,
                    center = TRUE,
                    shape = "QC",
                    color = "QC")
# Plot
pqn_pca + pqn_bp
```

![](Multivariate_Statistics_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

The data was transformed by the generalized logarithm (glog) method.

``` r
# Extract clean data
ppm_noflag <- drop_flagged(pqn_set)
# "SummarizedExperiment" package installation
#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("SummarizedExperiment")
library(SummarizedExperiment)
# Convert feature height table to SummarizedExperiment class
pmp_data <- SummarizedExperiment(assays = exprs(ppm_noflag),
                                 colData = ppm_noflag@phenoData@data)
# Package for generalized logarithmic transform
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("pmp")
library(pmp)
# Generalised logarithmic transform
glog_exprs <- glog_transformation(df = pmp_data@assays@data@listData[[1]],
                                  classes = pmp_data$QC,
                                  qc_label = "QC")
# Adding glog transformation to notame MetaboSet
glog_set <- ppm_noflag
exprs(glog_set) <- glog_exprs
```

Finally, the data is ready to be exported and proceed with the
statistical analysis.

``` r
#save(glog_set, file = paste0(ppath, "Result/notame_Result/notame_out.RData"))
```

# Univariate statistic

## Volcano plot

The volcano plot was implemented to support the multivariate statistic
by analyzing individual features between levels of light factor (light
and shade). The generalized logarithm (glog) transformed data, and the
probabilistic quotient normalization (PQN) data were used for volcano
plot analysis. The glog-transformed data was used for the two-sample
t-test and Welch’s t-test calculation. The two-sample t-test (for
features with equal variances) and Welch’s t-test (for features with
unequal variances) methods were used for the “y” axis of the volcano
plot, and fold change (calculated using PQN normalized data) was used
for the “x” axis. We plotted two volcanoes, one with SubQC data and the
other with individual data.

The first step was data extraction.

### SubQC volcano plot

SubQC data extraction for t-test calculation.

``` r
# Removal of sample group from the dataset
volc_glog <- glog_set[, glog_set$Analysis_group != "Sample"]
pData(volc_glog) <- droplevels(pData(volc_glog))
# Drop QC
volc_glog <- drop_qcs(volc_glog)
# Extracting SubQC light factor
volc_glog <- volc_glog[, volc_glog$Light_factor != "Other factor"]
pData(volc_glog) <- droplevels(pData(volc_glog))
```

SubQC data extraction for fold change calculation.

``` r
# Removal of sample group from the dataset
volc_pqn <- ppm_noflag[, ppm_noflag$Analysis_group != "Sample"]
pData(volc_pqn) <- droplevels(pData(volc_pqn))
# Drop QC
volc_pqn <- drop_qcs(volc_pqn)
# Extracting SubQC light factor
volc_pqn <- volc_pqn[, volc_pqn$Light_factor != "Other factor"]
pData(volc_pqn) <- droplevels(pData(volc_pqn))
```

Before the t-test calculation, the features homoscedasticity was tested
using Levene’s test. For features with equal variance, the Student’s
t-test was implemented, and for features with unequal variance, Welch’s
t-test was implemented.

``` r
# Performing homoscedasticity test
volc_th <-
  perform_homoscedasticity_tests(volc_glog,
                                 formula_char = "Feature ~ Light_factor")
```

    ## INFO [2025-01-31 19:31:03] Starting homoscedasticity tests.
    ## INFO [2025-01-31 19:31:20] Homoscedasticity tests performed.

According to Levene’s test, all features showed an equal variance
(p-value \> 0.05). Thus, the t-test and fold change will be calculated.

``` r
# Library to left_join use
library(dplyr)
# Fold change between shade and light levels
volc_fc <- fold_change(volc_pqn, group = "Light_factor")
```

    ## INFO [2025-01-31 19:31:20] Starting to compute fold changes.
    ## INFO [2025-01-31 19:31:21] Fold changes computed.

``` r
# two-sample t-test performing
volc_t <- perform_t_test(volc_glog,
                         formula_char = "Feature ~ Light_factor",
                         var.equal = TRUE)
```

    ## INFO [2025-01-31 19:31:21] Starting t-tests for Light & Shade
    ## INFO [2025-01-31 19:31:27] t-tests performed.

``` r
# Adding the fold change to the t-test data
volc_data <- left_join(volc_t, volc_fc)
# Log-transform for visualization
volc_data$logP <- -log10(volc_data$Light_vs_Shade_t_test_P)
volc_data$log2_fc <- log2(volc_data$Shade_vs_Light_FC)
```

Volcano plot of SubQC.

``` r
# Determine point colors based on significance and fold change
volc_data <- volc_data %>%
  mutate(point_color = case_when(
    Light_vs_Shade_t_test_P < 0.05 & log2_fc < -1 ~ "Down", # significantly down
    Light_vs_Shade_t_test_P < 0.05 & log2_fc > 1 ~ "Up",    # significantly up
    log2_fc < 1 & log2_fc > -1 ~ "No fold change",          # fold change < 2
    Light_vs_Shade_t_test_P > 0.05 & log2_fc < -1 | 
      Light_vs_Shade_t_test_P > 0.05 & log2_fc > 1 ~ "Fold change")) # fold change > 2
# Extracting feature identified
metab_data <- glog_set[!is.na(glog_set@featureData@data$Metabolite),]
# Extracting metabolite table
meta_table <- metab_data@featureData@data
# Creating a new small table of the annotated compounds
# Keep metabolites with p-value < 0.05
volc_compouds <- subset(volc_data, Light_vs_Shade_t_test_P < 0.05 |
                             point_color == "Fold change")
volc_compouds <- left_join(meta_table, volc_compouds)
# Volcano plot
vc_plot <- ggplot(volc_data, aes(log2_fc, logP, color = point_color)) +
  geom_point(size = 2, alpha = 0.4) +
  scale_colour_manual(values = c("No fold change" = "#808080",
                                 "Down" = "#4F9D4EFF",
                                 "Up" = "#ff80ff",
                                 "Fold change" = "#F5BC5CFF")) +
  theme_classic() +
  #theme(legend.position = "none") +
  #geom_point(data = volc_compouds,
  #           aes(shape = meta_table$Identification_level,
  #               color = meta_table$Identification_level),
  #           size = 2) +
  labs(#shape = 'Identification level',
       color = 'Color key') +
  ggrepel::geom_label_repel(data = volc_compouds,
                            aes(label = meta_table$Metabolite),
                            color = "black",
                            box.padding = 0.37,
                            label.padding = 0.22,
                            label.r = 0.30,
                            cex = 2.5,
                            max.overlaps = 20,
                            min.segment.length = 0,
                            seed = 42,) +
  xlab(bquote(Log[2](Shade/Light))) + ylab(bquote(-Log[10](p-value))) +
  scale_y_continuous(limits = c(0,4), labels = function(i) 10^-i,
                     expand=c(0.003, 0.003)) +
  theme(legend.position = c(0.90, 0.947),
        legend.background = element_rect(fill = "white", color = "black")) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(fill= "transparent")) +
  #geom_vline(xintercept = 1, linetype = "longdash", colour="gray") +
  #geom_vline(xintercept = -1, linetype = "longdash", colour="gray") +
  geom_vline(xintercept = 0, linetype = "longdash", colour="gray") +
  geom_hline(yintercept = -log10(0.05), linetype = "longdash", colour="gray")
vc_plot
```

![](Multivariate_Statistics_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

``` r
# Save plot
#ggsave('Result/notame_Result/figure_2b.pdf',
#       width = 5, height = 4, device='pdf', dpi="print")
```

### Sample volcano plot

Sample data extraction for t-test calculation.

``` r
# Removal of sample group from the dataset
sp_volc_glog <- glog_set[, glog_set$Analysis_group != "Sub_Quality_Control"]
pData(sp_volc_glog) <- droplevels(pData(sp_volc_glog))
# Drop QC
sp_volc_glog <- drop_qcs(sp_volc_glog)
```

Sample data extraction for fold change calculation.

``` r
# Removal of sample group from the dataset
sp_volc_pqn <- ppm_noflag[, ppm_noflag$Analysis_group != "Sub_Quality_Control"]
pData(sp_volc_pqn) <- droplevels(pData(sp_volc_pqn))
# Drop QC
sp_volc_pqn <- drop_qcs(sp_volc_pqn)
```

Before the t-test calculation, the features homoscedasticity was tested
using Levene’s test. For features with equal variance, the Student’s
t-test was implemented, and for features with unequal variance, Welch’s
t-test was implemented.

``` r
# Performing homoscedasticity test
sp_volc_th <-
  perform_homoscedasticity_tests(sp_volc_glog,
                                 formula_char = "Feature ~ Light_factor")
```

    ## INFO [2025-01-31 19:31:30] Starting homoscedasticity tests.
    ## INFO [2025-01-31 19:31:46] Homoscedasticity tests performed.

``` r
# Adding homoscedasticity results to notame MetaboSet
sp_volc_glog <- join_fData(sp_volc_glog, sp_volc_th)
# Extracting features with equal variance
sp_volc_ttset <-
  sp_volc_glog[sp_volc_glog@featureData@data$Levene_P > 0.05,]
# Extracting features with unequal variance
sp_volc_wttset <-
  sp_volc_glog[sp_volc_glog@featureData@data$Levene_P < 0.05,]
```

Calculation of t-test, Welch’s t-test and fold change.

``` r
# Fold change between shade and light levels
sp_volc_fc <- fold_change(sp_volc_pqn, group = "Light_factor")
```

    ## INFO [2025-01-31 19:31:47] Starting to compute fold changes.
    ## INFO [2025-01-31 19:31:47] Fold changes computed.

``` r
# The two-sample t-test performing
sp_volc_tt <- perform_t_test(sp_volc_ttset,
                            formula_char = "Feature ~ Light_factor",
                            var.equal = TRUE)
```

    ## INFO [2025-01-31 19:31:47] Starting t-tests for Light & Shade
    ## INFO [2025-01-31 19:31:53] t-tests performed.

``` r
# Adding a tag for t-test results
sp_volc_tt <- sp_volc_tt %>%
  mutate(Statistic_test = case_when(
    Feature_ID != 0 ~ "Student's t-test"))
# The Welch's t-test performing
sp_volc_wtt <- perform_t_test(sp_volc_wttset,
                            formula_char = "Feature ~ Light_factor",
                            var.equal = TRUE)
```

    ## INFO [2025-01-31 19:31:53] Starting t-tests for Light & Shade
    ## INFO [2025-01-31 19:31:54] t-tests performed.

``` r
# Adding a tag for Welch's t-test results
sp_volc_wtt <- sp_volc_wtt %>%
  mutate(Statistic_test = case_when(
    Feature_ID != 0 ~ "Welch's t-test"))
# Merge Welch's t-test and Student's t-test
sp_volc_t <- rbind(sp_volc_tt, sp_volc_wtt)
# Adding the fold change to statistical results
sp_volc_data <- left_join(sp_volc_t, sp_volc_fc)
# Log-transform for visualization
sp_volc_data$logP <- -log10(sp_volc_data$Light_vs_Shade_t_test_P)
sp_volc_data$log2_fc <- log2(sp_volc_data$Shade_vs_Light_FC)
```

Volcano plot of Sample.

``` r
# Determine point colors based on significance and fold change
sp_volc_data <- sp_volc_data %>%
  mutate(point_color = case_when(
    Light_vs_Shade_t_test_P < 0.05 & log2_fc < -1 ~ "Down", # significantly down
    Light_vs_Shade_t_test_P < 0.05 & log2_fc > 1 ~ "Up",    # significantly up
    log2_fc < 1 & log2_fc > -1 ~ "No fold change",          # fold change < 2
    Light_vs_Shade_t_test_P > 0.05 & log2_fc < -1 | 
      Light_vs_Shade_t_test_P > 0.05 & log2_fc > 1 ~ "Fold change")) # fold change > 2
# Creating a new small table of the annotated compounds
# Keep metabolites with p-value < 0.05
sp_volc_compouds <- subset(sp_volc_data, Light_vs_Shade_t_test_P < 0.05 |
                             point_color == "Fold change")
sp_volc_compouds <- left_join(meta_table, sp_volc_compouds)
# Volcano plot
sp_vc_plot <- ggplot(sp_volc_data, aes(log2_fc, logP)) +
  #geom_point(size = 2, alpha = 0.4) +
  scale_colour_manual(values = c("No fold change" = "#808080",
                                 "Down" = "#4F9D4EFF",
                                 "Up" = "#ff80ff",
                                 "Fold change" = "#F5BC5CFF")) +
  theme_classic() +
  #theme(legend.position = "none") +
  geom_point(data = sp_volc_data,
             aes(shape = Statistic_test,
                 color = point_color),
             size = 2.5,
             alpha = 0.4) +
  labs(shape = 'Statistical test',
       color = 'Color key') +
  ggrepel::geom_label_repel(data = sp_volc_compouds,
                            aes(label = meta_table$Metabolite),
                            color = "black",
                            box.padding = 0.37,
                            label.padding = 0.22,
                            label.r = 0.30,
                            cex = 2.5,
                            max.overlaps = 20,
                            min.segment.length = 0,
                            seed = 42) +
  xlab(bquote(Log[2](Shade/Light))) + ylab(bquote(-Log[10](p-value))) +
  scale_y_continuous(limits = c(0,4), labels = function(i) 10^-i,
                     expand=c(0.003, 0.003)) +
  theme(legend.position = c(0.915, 0.85),
        legend.background = element_rect(fill = "white", color = "black")) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(fill= "transparent")) +
  #geom_vline(xintercept = 1, linetype = "longdash", colour="gray") +
  #geom_vline(xintercept = -1, linetype = "longdash", colour="gray") +
  geom_vline(xintercept = 0, linetype = "longdash", colour="gray") +
  geom_hline(yintercept = -log10(0.05), linetype = "longdash", colour="gray")
sp_vc_plot
```

![](Multivariate_Statistics_files/figure-gfm/unnamed-chunk-28-1.png)<!-- -->

``` r
# Save plot
#ggsave(filename = "Result/notame_Result/Figure_s5.pdf", plot = sp_vc_plot,
#       width = 6, height = 5.4, units = "in", dpi = 300, scale = 1.7)
#ggsave(filename = "Result/notame_Result/Figure_s5.png", plot = sp_vc_plot,
#       width = 6, height = 5.4, units = "in", dpi = 300, scale = 1.7)
```

## Tukey test

The tukey test was implemented to support the multivariate statistic by
analyzing individual identified metabolites between levels of age factor
(early, medium, and late). Only metabolites with data homoscedasticity
were considered for the ANOVA test, and metabolites with significance
(p-value \< 0.05) were used for the Tukey test. The generalized
logarithm (glog) transformed data was used for the Tukey test.

Previously, we tested with the auto-scaled data and found that seven
metabolites did not satisfy the data homoscedasticity according to the
Bartlett test. On the other hand, when we tested with the
glog-transformed data, only two metabolites did not satisfy the data
homoscedasticity. For the above reasons, we used log-transformed data
for the Tukey test.

We analyzed two batches of data, one of the SubQC data and the other of
the individual data.

### SubQC Tukey test

SubQC data extraction for Tukey test calculation.

``` r
# Drop QC
tk_no_qc <- drop_qcs(glog_set)
# Drop samples
tk_no_sample <- tk_no_qc[, tk_no_qc$Analysis_group != "Sample"]
pData(tk_no_sample) <- droplevels(pData(tk_no_sample))
# Extracting levels of age factor
tk_age <- tk_no_sample[, tk_no_sample$Age_factor != "Other factor"]
pData(tk_age) <- droplevels(pData(tk_age))
# Extracting identified metabolites
tk_age_set <- tk_age[!is.na(tk_age@featureData@data$Metabolite),]
```

Preparing data and metadata for Tukey test calculation.

``` r
# Extracting feature height table
tk_height <- exprs(tk_age_set)
# Extracting sample information
tk_pdata <- tk_age_set@phenoData@data
tk_pdata <- data.frame(Sample_ID = tk_pdata$Sample_ID,
                       Age_factor = tk_pdata$Age_factor)
# Extracting feature information
tk_fdata <- tk_age_set@featureData@data
tk_fdata <- data.frame(Feature_ID = tk_fdata$Feature_ID,
                       Metabolite = tk_fdata$Metabolite)
# Transposing feature height table
tk_height <- t(tk_height)
# Feature height table to dataframe
tk_height <- as.data.frame(tk_height)
# Convert the row names to sample ID
tk_height <- tk_height %>% 
  mutate(Sample_ID = rownames(tk_height))
# Adding age factor as variable
tk_height <- left_join(tk_pdata, tk_height)
```

    ## Joining with `by = join_by(Sample_ID)`

The data is ready to test the data homoscedasticity. The Levene’s test
was used to inspect the homogeneity of variance.

``` r
# Performing homoscedasticity test
tk_bar_res <-
  perform_homoscedasticity_tests(tk_age_set,
                                 formula_char = "Feature ~ Age_factor")
```

    ## INFO [2025-01-31 19:31:56] Starting homoscedasticity tests.
    ## INFO [2025-01-31 19:31:57] Homoscedasticity tests performed.

``` r
# Adding a description of the p-value
tk_bar_res <- tk_bar_res %>%
  mutate(bartlett_group = ifelse(Levene_P > 0.05,
                                 "equal variance",
                                 "unequal variance"))
```

Most metabolites showed equality of variances according to the Bartlett
test. The next step will be the ANOVA test.

``` r
# Create two empty vectors to add results
Feature_ID <- c()
anova_pValue <- c()
# Batch analysis
for(i in 3:80) {
  anova_model <- glm(tk_height[, i] ~ tk_height$Age_factor, data = tk_height)
  anova_res <- anova(anova_model, test = "LRT")
  anova_pValue <- c(anova_pValue, anova_res$`Pr(>Chi)`[2])
  Feature_ID <- c(Feature_ID, names(tk_height[i]))
  tk_anova_res <- data.frame(Feature_ID, anova_pValue)
}
# Adding a description of the p-value
tk_anova_res <- tk_anova_res %>%
  mutate(anova_group = ifelse(anova_pValue > 0.05,
                              "Not Significant",
                              "Significant"))
```

Consequently, the features with significance (p-value \< 0.05) according
to ANOVA were used for the Tukey test. Tukey test groups the different
variables (in this case, the levels of age factor: early, medium, and
late) using significance tests.

``` r
# agricolae package installation and library loadding
#install.packages("agricolae", repos = "https://cran.r-project.org")
library(agricolae)
# Create two empty vectors to add results
Feature_ID <- c()
tk_group_early <- c()
tk_group_late <- c()
tk_group_medium <- c()
# Batch analysis
for(i in 3:80) {
  tk_model <- aov(tk_height[, i] ~ tk_height$Age_factor, data = tk_height)
  tk_tes <- HSD.test(tk_model, "tk_height$Age_factor", group = TRUE)
  tk_group_early <- c(tk_group_early, tk_tes[["groups"]]["Early",2])
  tk_group_late <- c(tk_group_late, tk_tes[["groups"]]["Late",2])
  tk_group_medium <- c(tk_group_medium, tk_tes[["groups"]]["Medium",2])
  Feature_ID <- c(Feature_ID, names(tk_height[i]))
  tk_res <- data.frame(Feature_ID, tk_group_early, tk_group_late,
                       tk_group_medium)
  tk_res
}
```

Merge all results and add a flag in the Tukey result column: “-” for
metabolites with unequal variance according to the Bartlett test and
“ns” for metabolites not significant according to the ANOVA test.

``` r
# Merge all results obtained in the Tukey test
tk_result <- left_join(tk_bar_res, tk_anova_res)
```

    ## Joining with `by = join_by(Feature_ID)`

``` r
tk_result <- left_join(tk_result, tk_res)
```

    ## Joining with `by = join_by(Feature_ID)`

``` r
# Flag metabolites
for(i in 1:78) {
  # Flag "-" metabolites that has unequal variance
  if (tk_result$bartlett_group[i] == "unequal variance") {
        tk_result$tk_group_early[i] <- "-"
        tk_result$tk_group_medium[i] <- "-"
        tk_result$tk_group_late[i] <- "-"
  }
  # Flag "ns" metabolites not significant
  if (tk_result$anova_group[i] == "Not Significant" &
      tk_result$bartlett_group[i] == "equal variance") {
        tk_result$tk_group_early[i] <- "ns"
        tk_result$tk_group_medium[i] <- "ns"
        tk_result$tk_group_late[i] <- "ns"
}
 tk_result
}
```

Creating a Tukey result table to use in the heatmap.

``` r
# Tukey table to heatmap
hm_tk_table <- data.frame(Feature_ID = tk_result$Feature_ID,
                          Early = tk_result$tk_group_early,
                          Late = tk_result$tk_group_late,
                          Medium = tk_result$tk_group_medium)
```

Formatting a table to show Tukey’s result.

``` r
# The "gt" package installation
#install.packages("gt")
library(gt)
# Adding metabolite name
tk_table <- left_join(tk_fdata, hm_tk_table)
# Deleting the Feature_ID column
tk_table <- subset(tk_table, select = -c(Feature_ID) )
tk_table <- tk_table |>
  gt() |>
  fmt_markdown(
  columns = everything(),
  rows = everything(),
  md_engine = c("markdown", "commonmark")) |>
  tab_spanner(
    label = "Tukey test",
    columns = 2:4) |>
    tab_source_note(source_note = "Note: Different lowercase letter within the 
                  same row indicates significant differences based on Tukey's 
                  multiple comparisons (p-value < 0.05, n = 3). '-' represents 
                  features with unequal variance according to Levene's test. 
                  'ns' are the features that do not have statistical 
                  differences (p-value > 0.05) according to the ANOVA test.") |>
  tab_options(table.width = px(750)) |>
  as_raw_html()
```

Print Tukey results in a table.

``` r
tk_table
```

<div id="dcdlxhktxj" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
  
  <table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false" style="-webkit-font-smoothing: antialiased; -moz-osx-font-smoothing: grayscale; font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji'; display: table; border-collapse: collapse; line-height: normal; margin-left: auto; margin-right: auto; color: #333333; font-size: 16px; font-weight: normal; font-style: normal; background-color: #FFFFFF; width: 750px; border-top-style: solid; border-top-width: 2px; border-top-color: #A8A8A8; border-right-style: none; border-right-width: 2px; border-right-color: #D3D3D3; border-bottom-style: solid; border-bottom-width: 2px; border-bottom-color: #A8A8A8; border-left-style: none; border-left-width: 2px; border-left-color: #D3D3D3;" width="750" bgcolor="#FFFFFF">
  <thead style="border-style: none;">
    <tr class="gt_col_headings gt_spanner_row" style="border-style: none; border-top-style: solid; border-top-width: 2px; border-top-color: #D3D3D3; border-bottom-width: 2px; border-bottom-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; border-bottom-style: hidden;">
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="2" colspan="1" scope="col" id="Metabolite" style="border-style: none; color: #333333; background-color: #FFFFFF; font-size: 100%; font-weight: normal; text-transform: inherit; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: bottom; padding-top: 5px; padding-bottom: 6px; padding-left: 5px; padding-right: 5px; overflow-x: hidden; text-align: left;" bgcolor="#FFFFFF" valign="bottom" align="left">Metabolite</th>
      <th class="gt_center gt_columns_top_border gt_column_spanner_outer" rowspan="1" colspan="3" scope="colgroup" id="Tukey test" style="border-style: none; color: #333333; background-color: #FFFFFF; font-size: 100%; font-weight: normal; text-transform: inherit; padding-top: 0; padding-bottom: 0; padding-left: 4px; text-align: center; padding-right: 0;" bgcolor="#FFFFFF" align="center">
        <div class="gt_column_spanner" style="border-bottom-style: solid; border-bottom-width: 2px; border-bottom-color: #D3D3D3; vertical-align: bottom; padding-top: 5px; padding-bottom: 5px; overflow-x: hidden; display: inline-block; width: 100%;">Tukey test</div>
      </th>
    </tr>
    <tr class="gt_col_headings" style="border-style: none; border-top-style: solid; border-top-width: 2px; border-top-color: #D3D3D3; border-bottom-style: solid; border-bottom-width: 2px; border-bottom-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3;">
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="Early" style="border-style: none; color: #333333; background-color: #FFFFFF; font-size: 100%; font-weight: normal; text-transform: inherit; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: bottom; padding-top: 5px; padding-bottom: 6px; padding-left: 5px; padding-right: 5px; overflow-x: hidden; text-align: left;" bgcolor="#FFFFFF" valign="bottom" align="left">Early</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="Late" style="border-style: none; color: #333333; background-color: #FFFFFF; font-size: 100%; font-weight: normal; text-transform: inherit; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: bottom; padding-top: 5px; padding-bottom: 6px; padding-left: 5px; padding-right: 5px; overflow-x: hidden; text-align: left;" bgcolor="#FFFFFF" valign="bottom" align="left">Late</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="Medium" style="border-style: none; color: #333333; background-color: #FFFFFF; font-size: 100%; font-weight: normal; text-transform: inherit; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: bottom; padding-top: 5px; padding-bottom: 6px; padding-left: 5px; padding-right: 5px; overflow-x: hidden; text-align: left;" bgcolor="#FFFFFF" valign="bottom" align="left">Medium</th>
    </tr>
  </thead>
  <tbody class="gt_table_body" style="border-style: none; border-top-style: solid; border-top-width: 2px; border-top-color: #D3D3D3; border-bottom-style: solid; border-bottom-width: 2px; border-bottom-color: #D3D3D3;">
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">(E)-2,6-Dimethylocta-3,7-diene-2,6-diol</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">c</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">Coniferyl alcohol</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">Methyl cinnamate</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">Tagetone</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">1-Acetylbenzotriazole</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ns</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ns</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ns</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">1-Dodecanol</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">1-Hexadecanol</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">1-Hydroxy-3-(4-hydroxy-3-methoxyphenyl)acetone</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ns</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ns</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ns</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">1-Tetradecanol</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">2-Aminophenol</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">3-Formylindole</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ns</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ns</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ns</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">3-hydroxy-5,6-epoxy-.beta.-ionol</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">3-Hydroxy-5,6-epoxy-β-ionone</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">3-Hydroxy-β-damascone</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ab</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">3-Methylcatechol</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">4-Ethylcatechol</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ab</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">4-Methylcatechol</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">4-Vinylguaiacol|o-Acetyl-p-cresol</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ns</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ns</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ns</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">6-Hydroxy-4,4,7a-trimethyl-5,6,7,7a-tetrahydrobenzofuran-2(4H)-one</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ns</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ns</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ns</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">Acetosyringon</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ns</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ns</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ns</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">Benzeneacetaldehyde</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">c</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">Benzoic acid</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ns</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ns</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ns</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">Blumenol C</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">Butanoic acid, 2-methyl-, 2-methyl-2-propenyl ester</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">Caffeine</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ns</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ns</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ns</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">Canolol</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ab</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">Carbonic acid, octadecyl vinyl ester</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">Catechol</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">10-Heptadecenoic acid, methyl ester</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">Coumaran</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">Cyclohexenecarboxylic acid</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ns</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ns</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ns</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">Cyclotene</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">c</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">Decane</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">Dehydromevalonic lactone</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">Dihydro-3-oxo-β-ionol</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">Dihydroconiferyl alcohol</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ns</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ns</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ns</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">Dihydrosyringenin</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">Dill apiol</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ns</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ns</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ns</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">Diphenylacetylene</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">c</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">Docosane</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ns</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ns</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ns</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">Ethyl palmitate</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">Guaiacylacetone</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">Guajol</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">Heneicosane</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ns</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ns</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ns</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">Hexacosane</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">Hydrocinnamic acid</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">Hydroquinone| Cinnamaldehyde</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">Indole</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ab</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">Isobutyl phthalate</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ns</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ns</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ns</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">Levoglucosan</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ab</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">Methyl (2E)-3-(4-hydroxyphenyl)-2-propenoate</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">Methyl 1-methylpiperidine-2-carboxylate</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">Methyl 2-(4-hydroxy-3,5-dimethoxyphenyl)acetate</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">Methyl laurate</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">Methyl linolenate</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">Methyl oleate</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">Methyl palmitate</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">c</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">Methyl palmitelaidate</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">Methylethylmaleimide</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ns</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ns</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ns</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">Nonadecane</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ns</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ns</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ns</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">Octadecane</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">Palmitic acid</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">Pentacosane</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">Phenol</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">Phenylethyl Alcohol</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">c</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">Phytol</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">c</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">Proline, 2,2-dimethyl-5-oxo-</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">Senecioic acid</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ab</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">Squalene</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ns</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ns</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ns</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">Syringenin</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ns</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ns</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ns</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">Syringol</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ns</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ns</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ns</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">Tetracosane</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">Theobromine</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">Theophylline</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">Tiglic acid</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ns</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ns</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">ns</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">4-Propenylsyringol</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">Undecane</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td></tr>
    <tr style="border-style: none;"><td headers="Metabolite" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">Vanillin</span></td>
<td headers="Early" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">a</span></td>
<td headers="Late" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td>
<td headers="Medium" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left"><span class="gt_from_md">b</span></td></tr>
  </tbody>
  <tfoot class="gt_sourcenotes" style="border-style: none; color: #333333; background-color: #FFFFFF; border-bottom-style: none; border-bottom-width: 2px; border-bottom-color: #D3D3D3; border-left-style: none; border-left-width: 2px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 2px; border-right-color: #D3D3D3;" bgcolor="#FFFFFF">
    <tr style="border-style: none;">
      <td class="gt_sourcenote" colspan="4" style="border-style: none; font-size: 90%; padding-top: 4px; padding-bottom: 4px; padding-left: 5px; padding-right: 5px;">Note: Different lowercase letter within the 
                  same row indicates significant differences based on Tukey's 
                  multiple comparisons (p-value &lt; 0.05, n = 3). '-' represents 
                  features with unequal variance according to Levene's test. 
                  'ns' are the features that do not have statistical 
                  differences (p-value &gt; 0.05) according to the ANOVA test.</td>
    </tr>
  </tfoot>
  
</table>
</div>

### Sample Tukey test

Sample data extraction for Tukey test calculation.

``` r
# Drop factor pooled (subQC)
tk_no_subqc <- tk_no_qc[, tk_no_qc$Analysis_group != "Sub_Quality_Control"]
pData(tk_no_subqc) <- droplevels(pData(tk_no_subqc))
# Extracting identified metabolites
tk_no_subqc <- tk_no_subqc[!is.na(tk_no_subqc@featureData@data$Metabolite),]
```

Preparing data and metadata for Tukey test calculation.

``` r
# Extracting feature height table
tk_spl_height <- exprs(tk_no_subqc)
# Extracting sample information
tk_spl_pdata <- tk_no_subqc@phenoData@data
tk_spl_pdata <- data.frame(Sample_ID = tk_spl_pdata$Sample_ID,
                           Age_factor = tk_spl_pdata$Age_factor)
# Extracting feature information
tk_spl_fdata <- tk_no_subqc@featureData@data
tk_spl_fdata <- data.frame(Feature_ID = tk_spl_fdata$Feature_ID,
                           Metabolite = tk_spl_fdata$Metabolite)
# Transposing feature height table
tk_spl_height <- t(tk_spl_height)
# Feature height table to dataframe
tk_spl_height <- as.data.frame(tk_spl_height)
# Convert the row names to sample ID
tk_spl_height <- tk_spl_height %>% 
  mutate(Sample_ID = rownames(tk_spl_height))
# Adding age factor as variable
tk_spl_height <- left_join(tk_spl_pdata, tk_spl_height)
```

    ## Joining with `by = join_by(Sample_ID)`

The data is ready to test the data homoscedasticity. The Levene’s test
was used to inspect the homogeneity of variance.

``` r
# Performing homoscedasticity test
tk_spl_bar_res <-
  perform_homoscedasticity_tests(tk_no_subqc,
                                 formula_char = "Feature ~ Age_factor")
```

    ## INFO [2025-01-31 19:32:02] Starting homoscedasticity tests.
    ## INFO [2025-01-31 19:32:03] Homoscedasticity tests performed.

``` r
# Adding a description of the p-value
tk_spl_bar_res <- tk_spl_bar_res %>%
  mutate(bartlett_group = ifelse(Levene_P > 0.05,
                                 "equal variance",
                                 "unequal variance"))
```

Some metabolites showed equality of variances according to the Bartlett
test. The next step will be the ANOVA test.

``` r
# Create two empty vectors to add results
Feature_ID <- c()
anova_pValue <- c()
# Batch analysis
for(i in 3:80) {
  anova_model <- glm(tk_spl_height[, i] ~ tk_spl_height$Age_factor,
                     data = tk_spl_height)
  anova_res <- anova(anova_model, test = "LRT")
  anova_pValue <- c(anova_pValue, anova_res$`Pr(>Chi)`[2])
  Feature_ID <- c(Feature_ID, names(tk_spl_height[i]))
  tk_spl_anova_res <- data.frame(Feature_ID, anova_pValue)
}
# Adding a description of the p-value
tk_spl_anova_res <- tk_spl_anova_res %>%
  mutate(anova_group = ifelse(anova_pValue > 0.05,
                              "Not Significant",
                              "Significant"))
```

Consequently, the features with significance (p-value \< 0.05) according
to ANOVA were used for the Tukey test. Tukey test groups the different
variables (in this case, the levels of age factor: early, medium, and
late) using significance tests.

``` r
# Create two empty vectors to add results
Feature_ID <- c()
tk_group_early <- c()
tk_group_late <- c()
tk_group_medium <- c()
# Batch analysis
for(i in 3:80) {
  tk_model <- aov(tk_spl_height[, i] ~ tk_spl_height$Age_factor,
                  data = tk_spl_height)
  tk_tes <- HSD.test(tk_model, "tk_spl_height$Age_factor", group = TRUE)
  tk_group_early <- c(tk_group_early, tk_tes[["groups"]]["Early",2])
  tk_group_late <- c(tk_group_late, tk_tes[["groups"]]["Late",2])
  tk_group_medium <- c(tk_group_medium, tk_tes[["groups"]]["Medium",2])
  Feature_ID <- c(Feature_ID, names(tk_spl_height[i]))
  tk_spl_res <- data.frame(Feature_ID, tk_group_early, tk_group_late,
                       tk_group_medium)
  tk_spl_res
}
```

Merge all results and add a flag in the Tukey result column: “-” for
metabolites with unequal variance according to the Bartlett test and
“ns” for metabolites not significant according to the ANOVA test.

``` r
# Merge all results obtained in the Tukey test
tk_spl_result <- left_join(tk_spl_bar_res, tk_spl_anova_res)
```

    ## Joining with `by = join_by(Feature_ID)`

``` r
tk_spl_result <- left_join(tk_spl_result, tk_spl_res)
```

    ## Joining with `by = join_by(Feature_ID)`

``` r
# Flag "-" metabolites that has unequal variance
for(i in 1:78) {
  if (tk_spl_result$bartlett_group[i] == "unequal variance") {
        tk_spl_result$tk_group_early[i] <- "-"
        tk_spl_result$tk_group_medium[i] <- "-"
        tk_spl_result$tk_group_late[i] <- "-"
  }
  if (tk_spl_result$anova_group[i] == "Not Significant" &
      tk_spl_result$bartlett_group[i] == "equal variance") {
        tk_spl_result$tk_group_early[i] <- "ns"
        tk_spl_result$tk_group_medium[i] <- "ns"
        tk_spl_result$tk_group_late[i] <- "ns"
}
 tk_spl_result
}
```

Creating a Tukey result table to use in the heatmap.

``` r
# Tukey table to heatmap
hm_tk_table_spl <- data.frame(Feature_ID = tk_spl_result$Feature_ID,
                          Early = tk_spl_result$tk_group_early,
                          Late = tk_spl_result$tk_group_late,
                          Medium = tk_spl_result$tk_group_medium)
```

Print Tukey results in a table.

``` r
# Adding metabolite name
tk_table_spl <- left_join(tk_fdata, hm_tk_table_spl)
# Deleting the Feature_ID column
tk_table_spl <- subset(tk_table_spl, select = -c(Feature_ID) )
tk_table_spl <- tk_table_spl %>%
  gt()%>%
  tab_spanner(
    label = "Tukey test",
    columns = 2:4)%>%
    tab_source_note(source_note = "Note: Different lowercase letter within the 
                  same row indicates significant differences based on Tukey's 
                  multiple comparisons (p-value < 0.05, n = 18). '-' represents 
                  features with unequal variance according to Levene's test. 
                  'ns' are the features that do not have statistical 
                  differences (p-value > 0.05) according to the ANOVA test.")
tk_table_spl
```

<div id="mrvfytwduv" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#mrvfytwduv table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}

#mrvfytwduv thead, #mrvfytwduv tbody, #mrvfytwduv tfoot, #mrvfytwduv tr, #mrvfytwduv td, #mrvfytwduv th {
  border-style: none;
}

#mrvfytwduv p {
  margin: 0;
  padding: 0;
}

#mrvfytwduv .gt_table {
  display: table;
  border-collapse: collapse;
  line-height: normal;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}

#mrvfytwduv .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}

#mrvfytwduv .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#mrvfytwduv .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 3px;
  padding-bottom: 5px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#mrvfytwduv .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#mrvfytwduv .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#mrvfytwduv .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#mrvfytwduv .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}

#mrvfytwduv .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}

#mrvfytwduv .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#mrvfytwduv .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#mrvfytwduv .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}

#mrvfytwduv .gt_spanner_row {
  border-bottom-style: hidden;
}

#mrvfytwduv .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}

#mrvfytwduv .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}

#mrvfytwduv .gt_from_md > :first-child {
  margin-top: 0;
}

#mrvfytwduv .gt_from_md > :last-child {
  margin-bottom: 0;
}

#mrvfytwduv .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}

#mrvfytwduv .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}

#mrvfytwduv .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}

#mrvfytwduv .gt_row_group_first td {
  border-top-width: 2px;
}

#mrvfytwduv .gt_row_group_first th {
  border-top-width: 2px;
}

#mrvfytwduv .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#mrvfytwduv .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}

#mrvfytwduv .gt_first_summary_row.thick {
  border-top-width: 2px;
}

#mrvfytwduv .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#mrvfytwduv .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#mrvfytwduv .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#mrvfytwduv .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}

#mrvfytwduv .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#mrvfytwduv .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#mrvfytwduv .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#mrvfytwduv .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#mrvfytwduv .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#mrvfytwduv .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#mrvfytwduv .gt_left {
  text-align: left;
}

#mrvfytwduv .gt_center {
  text-align: center;
}

#mrvfytwduv .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#mrvfytwduv .gt_font_normal {
  font-weight: normal;
}

#mrvfytwduv .gt_font_bold {
  font-weight: bold;
}

#mrvfytwduv .gt_font_italic {
  font-style: italic;
}

#mrvfytwduv .gt_super {
  font-size: 65%;
}

#mrvfytwduv .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}

#mrvfytwduv .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}

#mrvfytwduv .gt_indent_1 {
  text-indent: 5px;
}

#mrvfytwduv .gt_indent_2 {
  text-indent: 10px;
}

#mrvfytwduv .gt_indent_3 {
  text-indent: 15px;
}

#mrvfytwduv .gt_indent_4 {
  text-indent: 20px;
}

#mrvfytwduv .gt_indent_5 {
  text-indent: 25px;
}

#mrvfytwduv .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}

#mrvfytwduv div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <thead>
    <tr class="gt_col_headings gt_spanner_row">
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="2" colspan="1" scope="col" id="Metabolite">Metabolite</th>
      <th class="gt_center gt_columns_top_border gt_column_spanner_outer" rowspan="1" colspan="3" scope="colgroup" id="Tukey test">
        <div class="gt_column_spanner">Tukey test</div>
      </th>
    </tr>
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="Early">Early</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="Late">Late</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="Medium">Medium</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td headers="Metabolite" class="gt_row gt_left">(E)-2,6-Dimethylocta-3,7-diene-2,6-diol</td>
<td headers="Early" class="gt_row gt_left">ns</td>
<td headers="Late" class="gt_row gt_left">ns</td>
<td headers="Medium" class="gt_row gt_left">ns</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">Coniferyl alcohol</td>
<td headers="Early" class="gt_row gt_left">b</td>
<td headers="Late" class="gt_row gt_left">a</td>
<td headers="Medium" class="gt_row gt_left">b</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">Methyl cinnamate</td>
<td headers="Early" class="gt_row gt_left">ns</td>
<td headers="Late" class="gt_row gt_left">ns</td>
<td headers="Medium" class="gt_row gt_left">ns</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">Tagetone</td>
<td headers="Early" class="gt_row gt_left">ns</td>
<td headers="Late" class="gt_row gt_left">ns</td>
<td headers="Medium" class="gt_row gt_left">ns</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">1-Acetylbenzotriazole</td>
<td headers="Early" class="gt_row gt_left">ns</td>
<td headers="Late" class="gt_row gt_left">ns</td>
<td headers="Medium" class="gt_row gt_left">ns</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">1-Dodecanol</td>
<td headers="Early" class="gt_row gt_left">ns</td>
<td headers="Late" class="gt_row gt_left">ns</td>
<td headers="Medium" class="gt_row gt_left">ns</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">1-Hexadecanol</td>
<td headers="Early" class="gt_row gt_left">ns</td>
<td headers="Late" class="gt_row gt_left">ns</td>
<td headers="Medium" class="gt_row gt_left">ns</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">1-Hydroxy-3-(4-hydroxy-3-methoxyphenyl)acetone</td>
<td headers="Early" class="gt_row gt_left">-</td>
<td headers="Late" class="gt_row gt_left">-</td>
<td headers="Medium" class="gt_row gt_left">-</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">1-Tetradecanol</td>
<td headers="Early" class="gt_row gt_left">ns</td>
<td headers="Late" class="gt_row gt_left">ns</td>
<td headers="Medium" class="gt_row gt_left">ns</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">2-Aminophenol</td>
<td headers="Early" class="gt_row gt_left">-</td>
<td headers="Late" class="gt_row gt_left">-</td>
<td headers="Medium" class="gt_row gt_left">-</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">3-Formylindole</td>
<td headers="Early" class="gt_row gt_left">-</td>
<td headers="Late" class="gt_row gt_left">-</td>
<td headers="Medium" class="gt_row gt_left">-</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">3-hydroxy-5,6-epoxy-.beta.-ionol</td>
<td headers="Early" class="gt_row gt_left">-</td>
<td headers="Late" class="gt_row gt_left">-</td>
<td headers="Medium" class="gt_row gt_left">-</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">3-Hydroxy-5,6-epoxy-β-ionone</td>
<td headers="Early" class="gt_row gt_left">ns</td>
<td headers="Late" class="gt_row gt_left">ns</td>
<td headers="Medium" class="gt_row gt_left">ns</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">3-Hydroxy-β-damascone</td>
<td headers="Early" class="gt_row gt_left">a</td>
<td headers="Late" class="gt_row gt_left">a</td>
<td headers="Medium" class="gt_row gt_left">a</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">3-Methylcatechol</td>
<td headers="Early" class="gt_row gt_left">-</td>
<td headers="Late" class="gt_row gt_left">-</td>
<td headers="Medium" class="gt_row gt_left">-</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">4-Ethylcatechol</td>
<td headers="Early" class="gt_row gt_left">ns</td>
<td headers="Late" class="gt_row gt_left">ns</td>
<td headers="Medium" class="gt_row gt_left">ns</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">4-Methylcatechol</td>
<td headers="Early" class="gt_row gt_left">-</td>
<td headers="Late" class="gt_row gt_left">-</td>
<td headers="Medium" class="gt_row gt_left">-</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">4-Vinylguaiacol|o-Acetyl-p-cresol</td>
<td headers="Early" class="gt_row gt_left">-</td>
<td headers="Late" class="gt_row gt_left">-</td>
<td headers="Medium" class="gt_row gt_left">-</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">6-Hydroxy-4,4,7a-trimethyl-5,6,7,7a-tetrahydrobenzofuran-2(4H)-one</td>
<td headers="Early" class="gt_row gt_left">-</td>
<td headers="Late" class="gt_row gt_left">-</td>
<td headers="Medium" class="gt_row gt_left">-</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">Acetosyringon</td>
<td headers="Early" class="gt_row gt_left">ns</td>
<td headers="Late" class="gt_row gt_left">ns</td>
<td headers="Medium" class="gt_row gt_left">ns</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">Benzeneacetaldehyde</td>
<td headers="Early" class="gt_row gt_left">ns</td>
<td headers="Late" class="gt_row gt_left">ns</td>
<td headers="Medium" class="gt_row gt_left">ns</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">Benzoic acid</td>
<td headers="Early" class="gt_row gt_left">ns</td>
<td headers="Late" class="gt_row gt_left">ns</td>
<td headers="Medium" class="gt_row gt_left">ns</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">Blumenol C</td>
<td headers="Early" class="gt_row gt_left">-</td>
<td headers="Late" class="gt_row gt_left">-</td>
<td headers="Medium" class="gt_row gt_left">-</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">Butanoic acid, 2-methyl-, 2-methyl-2-propenyl ester</td>
<td headers="Early" class="gt_row gt_left">ns</td>
<td headers="Late" class="gt_row gt_left">ns</td>
<td headers="Medium" class="gt_row gt_left">ns</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">Caffeine</td>
<td headers="Early" class="gt_row gt_left">ns</td>
<td headers="Late" class="gt_row gt_left">ns</td>
<td headers="Medium" class="gt_row gt_left">ns</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">Canolol</td>
<td headers="Early" class="gt_row gt_left">ns</td>
<td headers="Late" class="gt_row gt_left">ns</td>
<td headers="Medium" class="gt_row gt_left">ns</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">Carbonic acid, octadecyl vinyl ester</td>
<td headers="Early" class="gt_row gt_left">ns</td>
<td headers="Late" class="gt_row gt_left">ns</td>
<td headers="Medium" class="gt_row gt_left">ns</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">Catechol</td>
<td headers="Early" class="gt_row gt_left">-</td>
<td headers="Late" class="gt_row gt_left">-</td>
<td headers="Medium" class="gt_row gt_left">-</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">10-Heptadecenoic acid, methyl ester</td>
<td headers="Early" class="gt_row gt_left">ns</td>
<td headers="Late" class="gt_row gt_left">ns</td>
<td headers="Medium" class="gt_row gt_left">ns</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">Coumaran</td>
<td headers="Early" class="gt_row gt_left">a</td>
<td headers="Late" class="gt_row gt_left">b</td>
<td headers="Medium" class="gt_row gt_left">b</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">Cyclohexenecarboxylic acid</td>
<td headers="Early" class="gt_row gt_left">ns</td>
<td headers="Late" class="gt_row gt_left">ns</td>
<td headers="Medium" class="gt_row gt_left">ns</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">Cyclotene</td>
<td headers="Early" class="gt_row gt_left">b</td>
<td headers="Late" class="gt_row gt_left">a</td>
<td headers="Medium" class="gt_row gt_left">b</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">Decane</td>
<td headers="Early" class="gt_row gt_left">-</td>
<td headers="Late" class="gt_row gt_left">-</td>
<td headers="Medium" class="gt_row gt_left">-</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">Dehydromevalonic lactone</td>
<td headers="Early" class="gt_row gt_left">ns</td>
<td headers="Late" class="gt_row gt_left">ns</td>
<td headers="Medium" class="gt_row gt_left">ns</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">Dihydro-3-oxo-β-ionol</td>
<td headers="Early" class="gt_row gt_left">-</td>
<td headers="Late" class="gt_row gt_left">-</td>
<td headers="Medium" class="gt_row gt_left">-</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">Dihydroconiferyl alcohol</td>
<td headers="Early" class="gt_row gt_left">ns</td>
<td headers="Late" class="gt_row gt_left">ns</td>
<td headers="Medium" class="gt_row gt_left">ns</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">Dihydrosyringenin</td>
<td headers="Early" class="gt_row gt_left">-</td>
<td headers="Late" class="gt_row gt_left">-</td>
<td headers="Medium" class="gt_row gt_left">-</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">Dill apiol</td>
<td headers="Early" class="gt_row gt_left">ns</td>
<td headers="Late" class="gt_row gt_left">ns</td>
<td headers="Medium" class="gt_row gt_left">ns</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">Diphenylacetylene</td>
<td headers="Early" class="gt_row gt_left">ns</td>
<td headers="Late" class="gt_row gt_left">ns</td>
<td headers="Medium" class="gt_row gt_left">ns</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">Docosane</td>
<td headers="Early" class="gt_row gt_left">ns</td>
<td headers="Late" class="gt_row gt_left">ns</td>
<td headers="Medium" class="gt_row gt_left">ns</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">Ethyl palmitate</td>
<td headers="Early" class="gt_row gt_left">-</td>
<td headers="Late" class="gt_row gt_left">-</td>
<td headers="Medium" class="gt_row gt_left">-</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">Guaiacylacetone</td>
<td headers="Early" class="gt_row gt_left">ns</td>
<td headers="Late" class="gt_row gt_left">ns</td>
<td headers="Medium" class="gt_row gt_left">ns</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">Guajol</td>
<td headers="Early" class="gt_row gt_left">ns</td>
<td headers="Late" class="gt_row gt_left">ns</td>
<td headers="Medium" class="gt_row gt_left">ns</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">Heneicosane</td>
<td headers="Early" class="gt_row gt_left">ns</td>
<td headers="Late" class="gt_row gt_left">ns</td>
<td headers="Medium" class="gt_row gt_left">ns</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">Hexacosane</td>
<td headers="Early" class="gt_row gt_left">ns</td>
<td headers="Late" class="gt_row gt_left">ns</td>
<td headers="Medium" class="gt_row gt_left">ns</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">Hydrocinnamic acid</td>
<td headers="Early" class="gt_row gt_left">ns</td>
<td headers="Late" class="gt_row gt_left">ns</td>
<td headers="Medium" class="gt_row gt_left">ns</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">Hydroquinone| Cinnamaldehyde</td>
<td headers="Early" class="gt_row gt_left">a</td>
<td headers="Late" class="gt_row gt_left">a</td>
<td headers="Medium" class="gt_row gt_left">a</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">Indole</td>
<td headers="Early" class="gt_row gt_left">ns</td>
<td headers="Late" class="gt_row gt_left">ns</td>
<td headers="Medium" class="gt_row gt_left">ns</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">Isobutyl phthalate</td>
<td headers="Early" class="gt_row gt_left">ns</td>
<td headers="Late" class="gt_row gt_left">ns</td>
<td headers="Medium" class="gt_row gt_left">ns</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">Levoglucosan</td>
<td headers="Early" class="gt_row gt_left">-</td>
<td headers="Late" class="gt_row gt_left">-</td>
<td headers="Medium" class="gt_row gt_left">-</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">Methyl (2E)-3-(4-hydroxyphenyl)-2-propenoate</td>
<td headers="Early" class="gt_row gt_left">ns</td>
<td headers="Late" class="gt_row gt_left">ns</td>
<td headers="Medium" class="gt_row gt_left">ns</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">Methyl 1-methylpiperidine-2-carboxylate</td>
<td headers="Early" class="gt_row gt_left">-</td>
<td headers="Late" class="gt_row gt_left">-</td>
<td headers="Medium" class="gt_row gt_left">-</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">Methyl 2-(4-hydroxy-3,5-dimethoxyphenyl)acetate</td>
<td headers="Early" class="gt_row gt_left">ns</td>
<td headers="Late" class="gt_row gt_left">ns</td>
<td headers="Medium" class="gt_row gt_left">ns</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">Methyl laurate</td>
<td headers="Early" class="gt_row gt_left">ns</td>
<td headers="Late" class="gt_row gt_left">ns</td>
<td headers="Medium" class="gt_row gt_left">ns</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">Methyl linolenate</td>
<td headers="Early" class="gt_row gt_left">ns</td>
<td headers="Late" class="gt_row gt_left">ns</td>
<td headers="Medium" class="gt_row gt_left">ns</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">Methyl oleate</td>
<td headers="Early" class="gt_row gt_left">ns</td>
<td headers="Late" class="gt_row gt_left">ns</td>
<td headers="Medium" class="gt_row gt_left">ns</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">Methyl palmitate</td>
<td headers="Early" class="gt_row gt_left">-</td>
<td headers="Late" class="gt_row gt_left">-</td>
<td headers="Medium" class="gt_row gt_left">-</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">Methyl palmitelaidate</td>
<td headers="Early" class="gt_row gt_left">-</td>
<td headers="Late" class="gt_row gt_left">-</td>
<td headers="Medium" class="gt_row gt_left">-</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">Methylethylmaleimide</td>
<td headers="Early" class="gt_row gt_left">ns</td>
<td headers="Late" class="gt_row gt_left">ns</td>
<td headers="Medium" class="gt_row gt_left">ns</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">Nonadecane</td>
<td headers="Early" class="gt_row gt_left">ns</td>
<td headers="Late" class="gt_row gt_left">ns</td>
<td headers="Medium" class="gt_row gt_left">ns</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">Octadecane</td>
<td headers="Early" class="gt_row gt_left">ns</td>
<td headers="Late" class="gt_row gt_left">ns</td>
<td headers="Medium" class="gt_row gt_left">ns</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">Palmitic acid</td>
<td headers="Early" class="gt_row gt_left">b</td>
<td headers="Late" class="gt_row gt_left">b</td>
<td headers="Medium" class="gt_row gt_left">a</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">Pentacosane</td>
<td headers="Early" class="gt_row gt_left">ns</td>
<td headers="Late" class="gt_row gt_left">ns</td>
<td headers="Medium" class="gt_row gt_left">ns</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">Phenol</td>
<td headers="Early" class="gt_row gt_left">-</td>
<td headers="Late" class="gt_row gt_left">-</td>
<td headers="Medium" class="gt_row gt_left">-</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">Phenylethyl Alcohol</td>
<td headers="Early" class="gt_row gt_left">ns</td>
<td headers="Late" class="gt_row gt_left">ns</td>
<td headers="Medium" class="gt_row gt_left">ns</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">Phytol</td>
<td headers="Early" class="gt_row gt_left">ns</td>
<td headers="Late" class="gt_row gt_left">ns</td>
<td headers="Medium" class="gt_row gt_left">ns</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">Proline, 2,2-dimethyl-5-oxo-</td>
<td headers="Early" class="gt_row gt_left">-</td>
<td headers="Late" class="gt_row gt_left">-</td>
<td headers="Medium" class="gt_row gt_left">-</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">Senecioic acid</td>
<td headers="Early" class="gt_row gt_left">-</td>
<td headers="Late" class="gt_row gt_left">-</td>
<td headers="Medium" class="gt_row gt_left">-</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">Squalene</td>
<td headers="Early" class="gt_row gt_left">ns</td>
<td headers="Late" class="gt_row gt_left">ns</td>
<td headers="Medium" class="gt_row gt_left">ns</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">Syringenin</td>
<td headers="Early" class="gt_row gt_left">ns</td>
<td headers="Late" class="gt_row gt_left">ns</td>
<td headers="Medium" class="gt_row gt_left">ns</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">Syringol</td>
<td headers="Early" class="gt_row gt_left">-</td>
<td headers="Late" class="gt_row gt_left">-</td>
<td headers="Medium" class="gt_row gt_left">-</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">Tetracosane</td>
<td headers="Early" class="gt_row gt_left">ns</td>
<td headers="Late" class="gt_row gt_left">ns</td>
<td headers="Medium" class="gt_row gt_left">ns</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">Theobromine</td>
<td headers="Early" class="gt_row gt_left">ns</td>
<td headers="Late" class="gt_row gt_left">ns</td>
<td headers="Medium" class="gt_row gt_left">ns</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">Theophylline</td>
<td headers="Early" class="gt_row gt_left">ns</td>
<td headers="Late" class="gt_row gt_left">ns</td>
<td headers="Medium" class="gt_row gt_left">ns</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">Tiglic acid</td>
<td headers="Early" class="gt_row gt_left">-</td>
<td headers="Late" class="gt_row gt_left">-</td>
<td headers="Medium" class="gt_row gt_left">-</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">4-Propenylsyringol</td>
<td headers="Early" class="gt_row gt_left">ns</td>
<td headers="Late" class="gt_row gt_left">ns</td>
<td headers="Medium" class="gt_row gt_left">ns</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">Undecane</td>
<td headers="Early" class="gt_row gt_left">ns</td>
<td headers="Late" class="gt_row gt_left">ns</td>
<td headers="Medium" class="gt_row gt_left">ns</td></tr>
    <tr><td headers="Metabolite" class="gt_row gt_left">Vanillin</td>
<td headers="Early" class="gt_row gt_left">ns</td>
<td headers="Late" class="gt_row gt_left">ns</td>
<td headers="Medium" class="gt_row gt_left">ns</td></tr>
  </tbody>
  <tfoot class="gt_sourcenotes">
    <tr>
      <td class="gt_sourcenote" colspan="4">Note: Different lowercase letter within the 
                  same row indicates significant differences based on Tukey's 
                  multiple comparisons (p-value &lt; 0.05, n = 18). '-' represents 
                  features with unequal variance according to Levene's test. 
                  'ns' are the features that do not have statistical 
                  differences (p-value &gt; 0.05) according to the ANOVA test.</td>
    </tr>
  </tfoot>
  
</table>
</div>

# Multivariate statistic

## Principal Component Analysis (PCA)

Extracting the data and metadata of notame MetaboSet for principal
component analysis (PCA) plot.

``` r
# Extracting feature height table
pca_peakheight <- exprs(glog_set)
# Extracting Phenotipic data
pca_phenodata <- glog_set@phenoData@data
```

Transposing feature table and preparing the PCA data.

``` r
# Transposing feature height table
transp_table  <- t(pca_peakheight)
# Centering features
ei_pca <- prcomp(transp_table, center = TRUE, scale. = FALSE)
```

### Light factor PCA

Sunlight exposures of the plant PCA plot.

``` r
# PCA scores
scores <- ei_pca$x %>%                   # Get PC coordinates
  data.frame %>%                         # Convert to data frames
  mutate(Sample_ID = rownames(.)) %>%    # Create a new column with the sample names
  left_join(pca_phenodata)               # Adding metadata
# Change the order of legend column of DataFrame
scores$Light_factor <- factor(scores$Light_factor,
                              levels = c("Light",
                                         "Shade",
                                         "Other factor",
                                         "QC"))
scores$Analysis_group <- factor(scores$Analysis_group,
                                levels = c("Sample",
                                           "Sub_Quality_Control",
                                           "QC"),
                                labels=c("Sample",
                                         "SubQC",
                                         "QC"))
# PCA plot
figure_1a <- ggplot(scores,
                     aes(PC1, PC2, shape = Analysis_group,
                         color = Light_factor)) +
  scale_color_manual(values=c("#F5BC5CFF",
                              "#32373AFF",
                              "#D3D9D9FF",
                              "#00BFC4")) +
  scale_shape_manual(values=c(1, 17, 15)) +
  geom_point(size = 3) +
  guides(x=guide_axis(title = "PC1 (18.35 %)"),
         y=guide_axis(title = "PC2 (15.20 %)"),
         shape = guide_legend(order = 2),
         color = guide_legend(order = 1)) +
  labs(shape = 'Analysis group', color= 'Light factor') +
  theme_classic() +
#  theme(legend.position = c(0.88, 0.80),
#        legend.background = element_rect(fill = "white", color = "black")) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(fill= "transparent")) +
  geom_vline(xintercept = 0, linetype = "longdash", colour="gray") +
  geom_hline(yintercept = 0, linetype = "longdash", colour="gray")
figure_1a
```

![](Multivariate_Statistics_files/figure-gfm/unnamed-chunk-48-1.png)<!-- -->

``` r
# Save plot
#ggsave('Result/notame_Result/Figure_1a.pdf',
#       width = 5, height = 4, device='pdf', dpi="print")
```

### Age factor PCA

Plant ages PCA plot.

``` r
# PCA scores
scores <- ei_pca$x %>%                   # Get PC coordinates
  data.frame %>%                         # Convert to data frames
  mutate(Sample_ID = rownames(.)) %>%    # Create a new column with the sample names
  left_join(pca_phenodata)               # Adding metadata
# Change the order of legend column of DataFrame
scores$Age_factor <- factor(scores$Age_factor,
                 levels = c("Early",
                            "Medium",
                            "Late",
                            "Other factor",
                            "QC"))
scores$Analysis_group <- factor(scores$Analysis_group,
                 levels = c("Sample",
                            "Sub_Quality_Control",
                            "QC"),
                 labels=c("Sample",
                          "SubQC",
                          "QC"))
# PCA plot
figure_1b <- ggplot(scores,
                     aes(PC1, PC2, shape = Analysis_group,
                         color = Age_factor)) +
  scale_color_manual(values=c("#E76BF3",
                              "#00B0F6",
                              "#F8766D",
                              "#D3D9D9FF",
                              "#00BFC4")) +
  scale_shape_manual(values=c(1, 17, 15)) + 
  geom_point(size = 3) +
  guides(x=guide_axis(title = "PC1 (18.35 %)"),
         y=guide_axis(title = "PC2 (15.20 %)"),
         shape = guide_legend(order = 2),
         color = guide_legend(order = 1)) +
  labs(shape = 'Analysis group', color= 'Age factor') +
  theme_classic() +
#  theme(legend.position = c(0.88, 0.77),
#        legend.background = element_rect(fill = "white", color = "black")) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(fill= "transparent")) +
  geom_vline(xintercept = 0, linetype = "longdash", colour="gray") +
  geom_hline(yintercept = 0, linetype = "longdash", colour="gray")
figure_1b
```

![](Multivariate_Statistics_files/figure-gfm/unnamed-chunk-49-1.png)<!-- -->

``` r
# Save plot
#ggsave('Result/notame_Result/Figure_1b.pdf',
#       width = 5, height = 4, device='pdf', dpi="print")
```

### Location factor PCA

Sample collected location PCA plot.

``` r
# PCA scores
scores <- ei_pca$x %>%                   # Get PC coordinates
  data.frame %>%                         # Convert to data frames
  mutate(Sample_ID = rownames(.)) %>%    # Create a new column with the sample names
  left_join(pca_phenodata)               # Adding metadata
# Change the order of legend column of DataFrame
scores$Location_factor <- factor(scores$Location_factor,
                 levels = c("Alto Pano",
                            "Alto Tena",
                            "Talag",
                            "Other factor",
                            "QC"))
scores$Analysis_group <- factor(scores$Analysis_group,
                 levels = c("Sample",
                            "Sub_Quality_Control",
                            "QC"),
                 labels=c("Sample",
                          "SubQC",
                          "QC"))
# PCA plot
figure_1c <- ggplot(scores,
                     aes(PC1, PC2, shape = Analysis_group,
                         color = Location_factor)) +
  scale_color_manual(values=c("#B24422FF",
                              "#4F9D4EFF",
                              "#4457A5FF",
                              "#D3D9D9FF",
                              "#00BFC4")) +
  scale_shape_manual(values=c(1, 17, 15)) + 
  geom_point(size = 3) +
  guides(x=guide_axis(title = "PC1 (18.35 %)"),
         y=guide_axis(title = "PC2 (15.20 %)")) +
  labs(shape = 'Analysis group', color= 'Location factor') +
  theme_classic() +
#  theme(legend.position = c(0.88, 0.77),
#        legend.background = element_rect(fill = "white", color = "black")) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(fill= "transparent")) +
  geom_vline(xintercept = 0, linetype = "longdash", colour="gray") +
  geom_hline(yintercept = 0, linetype = "longdash", colour="gray")
figure_1c
```

![](Multivariate_Statistics_files/figure-gfm/unnamed-chunk-50-1.png)<!-- -->

``` r
# Save plot
#ggsave('Result/notame_Result/Figure_1c.pdf',
#       width = 5, height = 4, device='pdf', dpi="print")
```

Plotting a loading PCA result.

``` r
loadings <- ei_pca$rotation %>%           # Extract loadings
  data.frame(Feature_ID = rownames(.))    # New column with feat name
```

Creating an artificial table with feature names and a compound column
(with identified metabolites).

``` r
# Creating a new small table with the annotated compounds
ei_compouds <- left_join(meta_table, loadings)
# Plotting results
figure_1d <- ggplot(loadings, aes(PC1, PC2)) +
  geom_point(alpha = 0.2) +
  theme_classic() + 
  geom_point(data = ei_compouds,
             aes(shape = meta_table$Identification_level,
                 color = meta_table$Identification_level),
             size = 2) +
  labs(shape = 'Identification level',
       color = 'Identification level') +
  ggrepel::geom_label_repel(data = ei_compouds,
                            aes(label = meta_table$Metabolite),
                            box.padding = 0.22,
                            label.padding = 0.22,
                            label.r = 0.3,
                            cex = 2.5,
                            max.overlaps = 7,
                            min.segment.length = 0) +
  guides(x=guide_axis(title = "PC1 (18.35 %)"),
         y=guide_axis(title = "PC2 (15.20 %)")) +
  theme(legend.position = c(0.945, 0.913),
        legend.background = element_rect(fill = "white", color = "black")) +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  geom_vline(xintercept = 0, linetype = "longdash", colour="gray") +
  geom_hline(yintercept = 0, linetype = "longdash", colour="gray") +
  ggsci::scale_color_aaas()
figure_1d
```

![](Multivariate_Statistics_files/figure-gfm/unnamed-chunk-52-1.png)<!-- -->

``` r
# Save plot
#ggsave('Result/notame_Result/figure_1d.pdf',
#       width = 20, height = 10.5, device='pdf', dpi="print")
```

### Alto Pano location PCA

Extracting the Alto Pano location data.

``` r
# Drop factor pooled (subQC)
no_subqc <- glog_set[, glog_set$Analysis_group != "Sub_Quality_Control"]
pData(no_subqc) <- droplevels(pData(no_subqc))
# Extracting "Alto Pano" (B) data
b_data <- no_subqc[, no_subqc$Location_factor != "Alto Tena"]
pData(b_data) <- droplevels(pData(b_data))
b_data <- b_data[, b_data$Location_factor != "Talag"]
pData(b_data) <- droplevels(pData(b_data))
# Extracting feature height table
b_peakheight <- exprs(b_data)
# Extracting Phenotipic data
b_phenodata <- b_data@phenoData@data
```

Transposing feature table and preparing the PCA data.

``` r
# Transposing feature height table
b_transp_table  <- t(b_peakheight)
# Centering and Scaling features
b_pca <- prcomp(b_transp_table, center = TRUE, scale. = FALSE)
```

Plotting the Alto Pano PCA result.

``` r
# PCA scores
b_scores <- b_pca$x %>%                  # Get PC coordinates
  data.frame %>%                         # Convert to data frames
  mutate(Sample_ID = rownames(.)) %>%    # Create a new column with the sample names
  left_join(b_phenodata)                 # Adding metadata
# Change the order of legend column of DataFrame
b_scores$Age_factor <- factor(b_scores$Age_factor,
                              levels = c("Early",
                                         "Medium",
                                         "Late",
                                         "QC"))
b_scores$Light_factor <- factor(b_scores$Light_factor,
                                levels = c("Light",
                                           "Shade",
                                           "QC"))
# PCA plot
figure_s1a <- ggplot(b_scores,
                     aes(PC1, PC2, shape = Age_factor,
                         color = Light_factor)) +
  geom_point(size = 3) +
  guides(x=guide_axis(title = "PC1 (32.40 %)"),
         y=guide_axis(title = "PC2 (21.22 %)")) +
  scale_color_manual(values=c("#F5BC5CFF",
                              "#32373AFF",
                              "#00BFC4")) +
  scale_shape_manual(values=c(16, 17, 15, 8)) + 
  labs(shape = 'Age factor',
       color='Light factor') +
  theme_classic() +
#  theme(legend.position = c(0.90, 0.65),
#        legend.background = element_rect(fill = "white", color = "black")) +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  geom_vline(xintercept = 0, linetype = "longdash", colour="gray") +
  geom_hline(yintercept = 0, linetype = "longdash", colour="gray")
figure_s1a
```

![](Multivariate_Statistics_files/figure-gfm/unnamed-chunk-55-1.png)<!-- -->

``` r
# Save plot
#ggsave('Result/notame_Result/figure_s1a.pdf',
#       width = 5, height = 4, device='pdf', dpi="print")
```

Plotting the Alto Pano loading PCA result.

``` r
b_loadings <- b_pca$rotation %>%          # Extract loadings
  data.frame(Feature_ID = rownames(.))    # New column with feat name
```

Creating an artificial table with feature names and a compound column
(with identified metabolites).

``` r
# Creating a new small table with the annotated compounds
b_compouds <- left_join(meta_table, b_loadings)
# Plotting results
figure_s1d <- ggplot(b_loadings, aes(PC1, PC2)) +
  geom_point(alpha = 0.2) +
  theme_classic() + 
  geom_point(data = b_compouds,
             aes(shape = meta_table$Identification_level,
                 color = meta_table$Identification_level),
             size = 2) +
  labs(shape = 'Identification level',
       color = 'Identification level') +
  ggrepel::geom_label_repel(data = b_compouds,
                            aes(label = meta_table$Metabolite),
                            box.padding = 0.37,
                            label.padding = 0.22,
                            label.r = 0.30,
                            cex = 2.5,
                            max.overlaps = 7) +
  guides(x=guide_axis(title = "PC1 (32.40 %)"),
         y=guide_axis(title = "PC2 (21.22 %)")) +
  theme(legend.position = c(0.945, 0.913),
        legend.background = element_rect(fill = "white", color = "black")) +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  geom_vline(xintercept = 0, linetype = "longdash", colour="gray") +
  geom_hline(yintercept = 0, linetype = "longdash", colour="gray") +
  ggsci::scale_color_aaas()
figure_s1d
```

![](Multivariate_Statistics_files/figure-gfm/unnamed-chunk-57-1.png)<!-- -->

``` r
# Save plot
#ggsave('Result/notame_Result/figure_s1d.pdf',
#       width = 5, height = 4, device='pdf', dpi="print")
```

### Alto Tena location PCA

Extracting the Alto Tena location data.

``` r
# Extracting "Alto Tena" (C) data
c_data <- no_subqc[, no_subqc$Location_factor != "Alto Pano"]
pData(c_data) <- droplevels(pData(c_data))
c_data <- c_data[, c_data$Location_factor != "Talag"]
pData(c_data) <- droplevels(pData(c_data))
# Extracting feature height table
c_peakheight <- exprs(c_data )
# Extracting Phenotipic data
c_phenodata <- c_data@phenoData@data
```

Transposing feature table and preparing the PCA data.

``` r
# Transposing feature height table
c_transptable  <- t(c_peakheight)
# Centering and Scaling features
c_pca <- prcomp(c_transptable, center = TRUE, scale. = FALSE)
```

Plotting the Alto Tena PCA result.

``` r
# PCA scores
c_scores <- c_pca$x %>%                  # Get PC coordinates
  data.frame %>%                         # Convert to data frames
  mutate(Sample_ID = rownames(.)) %>%    # Create a new column with the sample names
  left_join(c_phenodata)                 # Adding metadata
# Change the order of legend column of DataFrame
c_scores$Age_factor <- factor(c_scores$Age_factor,
                              levels = c("Early",
                                         "Medium",
                                         "Late",
                                         "QC"))
c_scores$Light_factor <- factor(c_scores$Light_factor,
                                levels = c("Light",
                                           "Shade",
                                           "QC"))
# PCA plot
figure_s1b <- ggplot(c_scores,
                    aes(PC1, PC2, shape = Age_factor,
                        color = Light_factor)) +
  geom_point(size = 3) +
  guides(x=guide_axis(title = "PC1 (24.53 %)"),
         y=guide_axis(title = "PC2 (16.96 %)")) +
  scale_color_manual(values=c("#F5BC5CFF",
                              "#32373AFF",
                              "#00BFC4")) +
  scale_shape_manual(values=c(16, 17, 15, 8)) + 
  labs(shape = 'Age factor',
       color = 'Light factor') +
  theme_classic() +
#  theme(legend.position = c(0.08, 0.65),
#        legend.background = element_rect(fill = "white", color = "black")) +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  geom_vline(xintercept = 0, linetype = "longdash", colour="gray") +
  geom_hline(yintercept = 0, linetype = "longdash", colour="gray")
figure_s1b
```

![](Multivariate_Statistics_files/figure-gfm/unnamed-chunk-60-1.png)<!-- -->

``` r
# Save plot
#ggsave('Result/notame_Result/figure_s1b.pdf',
#       width = 5, height = 4, device='pdf', dpi="print")
```

Plotting the Alto Tena loading PCA result.

``` r
c_loadings <- c_pca$rotation %>%          # Extract loadings
  data.frame(Feature_ID = rownames(.))    # New column with feat name
```

Creating an artificial table with feature names and a compound column
(with identified metabolites).

``` r
# Creating a new small table of the annotated compounds
c_compouds <- left_join(meta_table, c_loadings)
# Plotting results
figure_s1e <- ggplot(c_loadings, aes(PC1, PC2)) +
  geom_point(alpha = 0.2) +
  theme_classic() + 
  geom_point(data = c_compouds,
             aes(shape = meta_table$Identification_level,
                 color = meta_table$Identification_level),
             size = 2) +
  labs(shape = 'Identification level',
       color = 'Identification level') +
  ggrepel::geom_label_repel(data = c_compouds,
                            aes(label = meta_table$Metabolite),
                            box.padding = 0.37,
                            label.padding = 0.22,
                            label.r = 0.30,
                            cex = 2.5,
                            max.overlaps = 7) +
  guides(x=guide_axis(title = "PC1 (24.53 %)"),
         y=guide_axis(title = "PC2 (16.96 %)")) +
  theme(legend.position = c(0.945, 0.913),
        legend.background = element_rect(fill = "white", color = "black")) +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  geom_vline(xintercept = 0, linetype = "longdash", colour="gray") +
  geom_hline(yintercept = 0, linetype = "longdash", colour="gray") +
  ggsci::scale_color_aaas()
figure_s1e
```

![](Multivariate_Statistics_files/figure-gfm/unnamed-chunk-62-1.png)<!-- -->

``` r
# Save plot
#ggsave('Result/notame_Result/figure_s1e.pdf',
#       width = 5, height = 4, device='pdf', dpi="print")
```

### Talag location PCA

Extracting the Talag location data.

``` r
# Extracting "Talag" (A) data
a_data <- no_subqc[, no_subqc$Location_factor != "Alto Pano"]
pData(a_data) <- droplevels(pData(a_data))
a_data <- a_data[, a_data$Location_factor != "Alto Tena"]
pData(a_data) <- droplevels(pData(a_data))
# Extracting feature height table
a_peakheight <- exprs(a_data )
# Extracting Phenotipic data
a_phenodata <- a_data@phenoData@data
```

Transposing feature table and preparing the PCA data.

``` r
# Transposing feature height table
a_transptable  <- t(a_peakheight)
# Centering and Scaling features
a_pca <- prcomp(a_transptable, center = TRUE, scale. = FALSE)
```

Plotting the Talag PCA result.

``` r
# PCA scores
a_scores <- a_pca$x %>%                 # Get PC coordinates
  data.frame %>%                        # Convert to data frames
  mutate(Sample_ID = rownames(.)) %>%   # Create a new column with the sample names
  left_join(a_phenodata)                # Adding metadata
# Change the order of legend column of DataFrame
a_scores$Age_factor <- factor(a_scores$Age_factor,
                              levels = c("Early",
                                         "Medium",
                                         "Late",
                                         "QC"))
a_scores$Light_factor <- factor(a_scores$Light_factor,
                               levels = c("Light",
                                          "Shade",
                                          "QC"))
# PCA plot
figure_s1c <- ggplot(a_scores,
                    aes(PC1, PC2, shape = Age_factor,
                        color = Light_factor)) +
  geom_point(size = 3) +
  guides(x=guide_axis(title = "PC1 (31.23 %)"),
         y=guide_axis(title = "PC2 (21.11 %)")) +
  scale_color_manual(values=c("#F5BC5CFF",
                              "#32373AFF",
                              "#00BFC4")) +
  scale_shape_manual(values=c(16, 17, 15, 8)) + 
  labs(shape = 'Age factor',
       color = 'Light factor') +
  theme_classic() +
#  theme(legend.position = c(1.08, 0.65),
#        legend.background = element_rect(fill = "white", color = "black")) +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  geom_vline(xintercept = 0, linetype = "longdash", colour="gray") +
  geom_hline(yintercept = 0, linetype = "longdash", colour="gray")
figure_s1c
```

![](Multivariate_Statistics_files/figure-gfm/unnamed-chunk-65-1.png)<!-- -->

``` r
# Save plot
#ggsave('Result/notame_Result/figure_s1c.pdf',
#       width = 5, height = 4, device='pdf', dpi="print")
```

Plotting the Talag loading PCA result.

``` r
a_loadings <- a_pca$rotation %>%          # Extract loadings
  data.frame(Feature_ID = rownames(.))    # New column with feat name
```

Creating an artificial table with feature names and a compound column
(with identified metabolites).

``` r
# Creating a new small table of the annotated compounds
a_compouds <- left_join(meta_table, a_loadings)
# Plotting results
figure_s1f <- ggplot(a_loadings, aes(PC1, PC2)) +
  geom_point(alpha = 0.2) +
  theme_classic() + 
  geom_point(data = a_compouds,
             aes(shape = meta_table$Identification_level,
                 color = meta_table$Identification_level),
             size = 2) +
  labs(shape = 'Identification level',
       color = 'Identification level') +
  ggrepel::geom_label_repel(data = a_compouds,
                            aes(label = meta_table$Metabolite),
                            box.padding = 0.37,
                            label.padding = 0.22,
                            label.r = 0.30,
                            cex = 2.5,
                            max.overlaps = 7) +
  guides(x=guide_axis(title = "PC1 (31.23 %)"),
         y=guide_axis(title = "PC2 (21.11 %)")) +
  theme(legend.position = c(0.945, 0.913),
        legend.background = element_rect(fill = "white", color = "black")) +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  geom_vline(xintercept = 0, linetype = "longdash", colour="gray") +
  geom_hline(yintercept = 0, linetype = "longdash", colour="gray") +
  ggsci::scale_color_aaas()
figure_s1f
```

![](Multivariate_Statistics_files/figure-gfm/unnamed-chunk-67-1.png)<!-- -->

``` r
# Save plot
#ggsave('Result/notame_Result/figure_s1f.pdf',
#       width = 5, height = 4, device='pdf', dpi="print")
```

## Heatmap and HCA

Only the annotated features were used in the heatmap and the
hierarchical cluster analysis (HCA). Also, the PCA result showed more
clear separation of the levels in the SubQC than individual samples; for
this reason, the SubQC data was used in the heatmap analysis. For better
heatmap plot visualization, we scaled (by autoscaling method) the
glog-transformed data. The ComplexHeatmap package will be used for the
heatmap and the hierarchical cluster analysis (HCA).

ComplexHeatmap package and other dependency package installation.

``` r
# ComplexHeatmap package installation
#if (!requireNamespace("BiocManager", quietly=TRUE))
#    install.packages("BiocManager")
#BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)

# colorRamp2 package installation
#if (!requireNamespace("devtools", quietly = TRUE)) {
#  install.packages("devtools")
#}
#devtools::install_github("jokergoo/colorRamp2")
library(colorRamp2)

# cowplot package installation
#install.packages("cowplot")
library(cowplot)

# mdatools package installation
#install_github('svkucheryavski/mdatools')
library(mdatools)

# ClassyFire package installation
#remotes::install_github('aberHRML/classyfireR')
library(classyfireR)
```

The metabolites were classified using the
[ClassyFireR](https://doi.org/10.1186/s13321-016-0174-y) to add these
metabolite classifications to the heatmap plot.

``` r
# InChI key of the metabolites you want to classify
InChI_Keys <- c(Benzeneacetaldehyde = "DTUQWGWMVIHBKE-UHFFFAOYSA-N")
# Get classification
Classification_List <- purrr::map(InChI_Keys, get_classification)
```

    ## <U+2714> DTUQWGWMVIHBKE-UHFFFAOYSA-N

``` r
Classification_List
```

    ## $Benzeneacetaldehyde
    ## -- ClassyFire Object ------------------------------------- classyfireR v0.3.8 -- 
    ## Object Size: 10.9 Kb 
    ##  
    ## Info: 
    ## * InChIKey=DTUQWGWMVIHBKE-UHFFFAOYSA-N
    ##   
    ## * O=CCC1=CC=CC=C1
    ##   
    ## * Classification Version: 2.1
    ##   
    ## kingdom : Organic compounds
    ## \-superclass : Benzenoids
    ##   \-class : Benzene and substituted derivatives
    ##     \-subclass : Phenylacetaldehydes

### Heatmap and HCA of age

Scaling glog data by auto-Scaling method and age factor data extraction.

``` r
# Scaling by autoscaling method
hm_scl <- scale(t(glog_exprs), center = TRUE, scale = TRUE)
hm_scl <- t(hm_scl)
# Adding autoscaling data to notame MetaboSet
hm_scl_set <- ppm_noflag
exprs(hm_scl_set) <- hm_scl
# Drop QC
hm_no_qc <- drop_qcs(hm_scl_set)
# Drop samples
hm_no_sample <- hm_no_qc[, hm_no_qc$Analysis_group != "Sample"]
pData(hm_no_sample) <- droplevels(pData(hm_no_sample))
# Extracting levels of age factor
hm_age <- hm_no_sample[, hm_no_sample$Age_factor != "Other factor"]
pData(hm_age) <- droplevels(pData(hm_age))
# Extracting identified metabolites
hm_age_set <- hm_age[!is.na(hm_age@featureData@data$Metabolite),]
```

We improve heatmap visualization by adding the Tukey test as row
annotation. The Tukey test calculation code is available in this
notebook’s [Univariate statistic](#univariate-statistic) section.

``` r
# Adding tukey test result
hm_age_set <- join_fData(hm_age_set, hm_tk_table)
```

Row and top heatmap annotation.

``` r
set.seed(2571)
# Extracting feature height table
hm_height <- exprs(hm_age_set)
# Extracting sample information
hm_pdata <- hm_age_set@phenoData@data
# Extracting feature information
hm_fdata <- hm_age_set@featureData@data
# Adding row and column name
rownames(hm_height) <- hm_fdata$Metabolite
colnames(hm_height) <- hm_pdata$Group
# Metabolite superclass color
cols_metclass <- c("Organic oxygen compounds" = "#C8876EFF",
                   "Benzenoids" = "#2F533CFF",
                   "Phenylpropanoids and polyketides" = "#63764EFF",
                   "Organoheterocyclic compounds" = "#A87735FF",
                   "Lipids and lipid-like molecules" = "#F3C571FF",
                   "Organic acids and derivatives" = "#EEB043FF",
                   "Hydrocarbons" = "#E0D2B5FF")
# Add row anotation to HeatMap
hm_row_ann <- rowAnnotation(`Superclass` = hm_fdata$classyfireR_Superclass,
                            `Tukey test 2` = anno_text(hm_fdata$Late,
                                                       location = 0.5,
                                                       just = "center",
                                                       gp = gpar(fill = "#F8766D",
                                                                 col = "#666666"),
                                                       show_name = TRUE),
                            `Tukey test 1` = anno_text(hm_fdata$Medium,
                                                       location = 0.5,
                                                       just = "center",
                                                       gp = gpar(fill = "#00B0F6",
                                                                 col = "#666666"),
                                                       show_name = TRUE),
                            `Tukey test 0` = anno_text(hm_fdata$Early,
                                                       location = 0.5,
                                                       just = "center",
                                                       gp = gpar(fill = "#E76BF3",
                                                                 col = "#666666"),
                                                       show_name = TRUE),
                            col = list(`Superclass` = cols_metclass),
                            show_annotation_name = T,
                            show_legend = F)
# Factor levels color
cols_levels <- c("Shade" = "#32373AFF",
                 "Early" = "#E76BF3",
                 "Light" = "#F5BC5CFF",
                 "Late" = "#F8766D",
                 "Talag" = "#4457A5FF",
                 "Medium" = "#00B0F6",
                 "Alto Tena" = "#4F9D4EFF",
                 "Alto Pano" = "#B24422FF")
cols_light <- c("Light" = "#F5BC5CFF",
                "Shade" = "#32373AFF")
cols_age <- c("Early" = "#E76BF3",
              "Medium" = "#00B0F6",
              "Late" = "#F8766D")
cols_location <- c("Alto Pano" = "#B24422FF",
                   "Alto Tena" = "#4F9D4EFF",
                   "Talag" = "#4457A5FF")
# Add top anotation to HeatMap
top_info_ann <- HeatmapAnnotation(`Age factor` = hm_pdata$Group,
                                  col = list(`Age factor` = cols_levels),
                                  show_annotation_name = T,
                                  show_legend=F, 
                                  border = TRUE)
# Color scale
mycol <- colorRamp2(c(-2, 0, 2),
                    c("blue", "white", "red"))
# Heatmap matrix plotting
hm_plot <- Heatmap(hm_height,
                   col = mycol,
                   border_gp = grid::gpar(col = "black", lty = 0.05),
                   rect_gp = grid::gpar(col = "black", lwd = 0.75),
                   clustering_distance_columns = "euclidean",
                   clustering_method_columns = "complete",
                   top_annotation = top_info_ann,
                   right_annotation = hm_row_ann,
                   row_names_max_width = unit(10, "cm"),
                   show_heatmap_legend = FALSE,
                   row_km = 3, column_km = 2,
                   row_title = c("a", "b", "c"))
hm_plot
```

![](Multivariate_Statistics_files/figure-gfm/unnamed-chunk-72-1.png)<!-- -->

Adding legends to heatmap.

``` r
# Color scale legend
lgd1 <- Legend(col_fun = mycol,
               title = "glog abundance",
               direction = "horizontal" )
# Factor levels legend
lgd2 <- Legend(labels = c("Light (Pos)",
                          "Shade (Neg)"),
               legend_gp = gpar(fill = cols_light),
               title = "Light factor", ncol = 1)
lgd3 <- Legend(labels = c("Early (0)",
                          "Medium (1)",
                          "Late (2)"),
               legend_gp = gpar(fill = cols_age),
               title = "Age factor", ncol = 1)
lgd4 <- Legend(labels = c("Alto Pano (B)",
                          "Alto Tena (C)",
                          "Talag (A)"),
               legend_gp = gpar(fill = cols_location),
               title = "Location factor", ncol = 1)
# Metabolite class Legend
lgd5 <- Legend(labels = c(unique(hm_fdata$classyfireR_Superclass)),
               legend_gp = gpar(fill = cols_metclass),
               title = "Metabolite superclass", ncol = 3)
```

Heatmap plot.

``` r
set.seed(2571)
# Converting to ggplot
gg_heatmap <- grid.grabExpr(draw(hm_plot))
gg_heatmap <- ggpubr::as_ggplot(gg_heatmap)
# Legends
all_legends <- packLegend(lgd1,
                          #lgd2,
                          lgd3,
                          #lgd4,
                          lgd5,
                          direction = "horizontal")
gg_legend <- grid.grabExpr(draw(all_legends))
gg_legend_fn <- ggpubr::as_ggplot(gg_legend)
# Heatmap plot
gcms_hm <- plot_grid(gg_legend_fn,
                     gg_heatmap, ncol = 1,
                     rel_heights = c(0.045, 0.880))
gcms_hm
```

![](Multivariate_Statistics_files/figure-gfm/unnamed-chunk-74-1.png)<!-- -->

``` r
# Save heatmap plot
#ggsave(filename = "Result/notame_Result/Figure_2a.pdf", plot = gcms_hm,
#       width = 7, height = 10, units = "in", dpi = 300, scale = 1.7)
#ggsave(filename = "Result/notame_Result/Figure_2a.png", plot = gcms_hm,
#       width = 7, height = 10, units = "in", dpi = 300, scale = 1.7)
```

### Heatmap and HCA of sample

Using the scaled data of the age heatmap, we extract sample data of
identified metabolites.

``` r
# Drop subQC
hm_no_subqc <- hm_no_qc[, hm_no_qc$Analysis_group != "Sub_Quality_Control"]
pData(hm_no_subqc) <- droplevels(pData(hm_no_subqc))
# Extracting identified metabolites
hm_age_spl <- hm_no_subqc[!is.na(hm_no_subqc@featureData@data$Metabolite),]
```

We improve heatmap visualization by adding the Tukey test as row
annotation. The Tukey test calculation code is available in this
notebook’s [Univariate statistic](#univariate-statistic) section.

``` r
# Adding tukey test result
hm_age_spl <- join_fData(hm_age_spl, hm_tk_table_spl)
```

Row and top heatmap annotation.

``` r
set.seed(2571)
# Extracting feature height table
hm_height_spl <- exprs(hm_age_spl)
# Extracting sample information
hm_pdata_spl <- hm_age_spl@phenoData@data
# Extracting feature information
hm_fdata_spl <- hm_age_spl@featureData@data
# Adding row and column name
rownames(hm_height_spl) <- hm_fdata_spl$Metabolite
colnames(hm_height_spl) <- hm_pdata_spl$Group
# Add row anotation to HeatMap
hm_row_ann <- rowAnnotation(`Superclass` = hm_fdata_spl$classyfireR_Superclass,
                            `Tukey test 2` = anno_text(hm_fdata_spl$Late,
                                                       location = 0.5,
                                                       just = "center",
                                                       gp = gpar(fill = "#F8766D",
                                                                 col = "#666666"),
                                                       show_name = TRUE),
                            `Tukey test 1` = anno_text(hm_fdata_spl$Medium,
                                                       location = 0.5,
                                                       just = "center",
                                                       gp = gpar(fill = "#00B0F6",
                                                                 col = "#666666"),
                                                       show_name = TRUE),
                            `Tukey test 0` = anno_text(hm_fdata_spl$Early,
                                                       location = 0.5,
                                                       just = "center",
                                                       gp = gpar(fill = "#E76BF3",
                                                                 col = "#666666"),
                                                       show_name = TRUE),
                            col = list(`Superclass` = cols_metclass),
                            show_annotation_name = T,
                            show_legend = F)
# Factor levels color
cols_levels <- c("Shade" = "#32373AFF",
                 "Early" = "#E76BF3",
                 "Light" = "#F5BC5CFF",
                 "Late" = "#F8766D",
                 "Talag" = "#4457A5FF",
                 "Medium" = "#00B0F6",
                 "Alto Tena" = "#4F9D4EFF",
                 "Alto Pano" = "#B24422FF")
cols_light <- c("Light" = "#F5BC5CFF",
                "Shade" = "#32373AFF")
cols_age <- c("Early" = "#E76BF3",
              "Medium" = "#00B0F6",
              "Late" = "#F8766D")
cols_location <- c("Alto Pano" = "#B24422FF",
                   "Alto Tena" = "#4F9D4EFF",
                   "Talag" = "#4457A5FF")
# Add top anotation to HeatMap
top_info_ann <- HeatmapAnnotation(`Light factor` = hm_pdata_spl$Light_factor,
                                  `Age factor` = hm_pdata_spl$Age_factor,
                                  `Location factor` = hm_pdata_spl$Location_factor,
                                  col = list(`Light factor` = cols_light,
                                             `Age factor` = cols_age,
                                             `Location factor` = cols_location),
                                  show_annotation_name = T,
                                  show_legend=F, 
                                  border = TRUE)
# Color scale
mycol_spl <- colorRamp2(c(-4, 0, 4),
                    c("blue", "white", "red"))
# Heatmap matrix plotting
hm_plot_spl <- Heatmap(hm_height_spl,
                       col = mycol,
                       border_gp = grid::gpar(col = "black", lty = 0.05),
                       rect_gp = grid::gpar(col = "black", lwd = 0.75),
                       clustering_distance_columns = "euclidean",
                       clustering_method_columns = "complete",
                       top_annotation = top_info_ann,
                       right_annotation = hm_row_ann,
                       row_names_max_width = unit(10, "cm"),
                       show_heatmap_legend = FALSE,
                       row_km = 3, column_km = 2,
                       row_title = c("a", "b", "c"))
hm_plot_spl
```

![](Multivariate_Statistics_files/figure-gfm/unnamed-chunk-77-1.png)<!-- -->

Adding legends to heatmap.

``` r
# Color scale legend
lgd_hm_spl <- Legend(col_fun = mycol_spl,
                     title = "glog abundance",
                     direction = "horizontal" )
```

Heatmap plot.

``` r
set.seed(2571)
# Converting to ggplot
gg_heatmap_spl <- grid.grabExpr(draw(hm_plot_spl))
gg_heatmap_spl <- ggpubr::as_ggplot(gg_heatmap_spl)
# Legends
all_legends_spl <- packLegend(lgd_hm_spl,
                           lgd2,
                           lgd3, 
                           lgd4, 
                           lgd5,
                           direction = "horizontal")
gg_legend_spl <- grid.grabExpr(draw(all_legends_spl))
gg_legend_fn_spl <- ggpubr::as_ggplot(gg_legend_spl)
# Heatmap plot
gcms_hm_spl <- plot_grid(gg_legend_fn_spl,
                     gg_heatmap_spl, ncol = 1,
                     rel_heights = c(0.045, 0.880))
gcms_hm_spl
```

![](Multivariate_Statistics_files/figure-gfm/unnamed-chunk-79-1.png)<!-- -->

``` r
# Save heatmap plot
#ggsave(filename = "Result/notame_Result/Figure_s4.pdf", plot = gcms_hm_spl,
#       width = 8, height = 10, units = "in", dpi = 300, scale = 1.7)
#ggsave(filename = "Result/notame_Result/Figure_s4.png", plot = gcms_hm_spl,
#       width = 8, height = 10, units = "in", dpi = 300, scale = 1.7)
```

Finish a record.

``` r
finish_log()
```

    ## INFO [2025-01-31 19:32:29] Finished analysis. Fri Jan 31 19:32:29 2025
    ## Session info:
    ## 
    ## INFO [2025-01-31 19:32:29] R version 4.4.1 (2024-06-14 ucrt)
    ## INFO [2025-01-31 19:32:29] Platform: x86_64-w64-mingw32/x64
    ## INFO [2025-01-31 19:32:29] Running under: Windows 7 x64 (build 7601) Service Pack 1
    ## INFO [2025-01-31 19:32:29] 
    ## INFO [2025-01-31 19:32:29] Matrix products: default
    ## INFO [2025-01-31 19:32:29] 
    ## INFO [2025-01-31 19:32:29] 
    ## INFO [2025-01-31 19:32:29] locale:
    ## INFO [2025-01-31 19:32:29] [1] LC_COLLATE=English_United States.1252 
    ## INFO [2025-01-31 19:32:29] [2] LC_CTYPE=English_United States.1252   
    ## INFO [2025-01-31 19:32:29] [3] LC_MONETARY=English_United States.1252
    ## INFO [2025-01-31 19:32:29] [4] LC_NUMERIC=C                          
    ## INFO [2025-01-31 19:32:29] [5] LC_TIME=English_United States.1252    
    ## INFO [2025-01-31 19:32:29] 
    ## INFO [2025-01-31 19:32:29] time zone: America/Guayaquil
    ## INFO [2025-01-31 19:32:29] tzcode source: internal
    ## INFO [2025-01-31 19:32:29] 
    ## INFO [2025-01-31 19:32:29] attached base packages:
    ## INFO [2025-01-31 19:32:29] [1] grid      stats4    stats     graphics  grDevices utils     datasets 
    ## INFO [2025-01-31 19:32:29] [8] methods   base     
    ## INFO [2025-01-31 19:32:29] 
    ## INFO [2025-01-31 19:32:29] other attached packages:
    ## INFO [2025-01-31 19:32:29]  [1] classyfireR_0.3.8           mdatools_0.14.2            
    ## INFO [2025-01-31 19:32:29]  [3] cowplot_1.1.3               colorRamp2_0.1.0           
    ## INFO [2025-01-31 19:32:29]  [5] ComplexHeatmap_2.20.0       gt_0.11.1                  
    ## INFO [2025-01-31 19:32:29]  [7] agricolae_1.3-7             dplyr_1.1.4                
    ## INFO [2025-01-31 19:32:29]  [9] pmp_1.14.1                  SummarizedExperiment_1.34.0
    ## INFO [2025-01-31 19:32:29] [11] GenomicRanges_1.56.0        GenomeInfoDb_1.40.1        
    ## INFO [2025-01-31 19:32:29] [13] IRanges_2.38.0              S4Vectors_0.42.0           
    ## INFO [2025-01-31 19:32:29] [15] MatrixGenerics_1.16.0       matrixStats_1.3.0          
    ## INFO [2025-01-31 19:32:29] [17] patchwork_1.2.0             notame_0.3.1               
    ## INFO [2025-01-31 19:32:29] [19] magrittr_2.0.3              ggplot2_3.5.1              
    ## INFO [2025-01-31 19:32:29] [21] futile.logger_1.4.3         Biobase_2.64.0             
    ## INFO [2025-01-31 19:32:29] [23] BiocGenerics_0.50.0        
    ## INFO [2025-01-31 19:32:29] 
    ## INFO [2025-01-31 19:32:29] loaded via a namespace (and not attached):
    ## INFO [2025-01-31 19:32:29]   [1] RColorBrewer_1.1-3      sys_3.4.2               rstudioapi_0.16.0      
    ## INFO [2025-01-31 19:32:29]   [4] jsonlite_1.8.8          shape_1.4.6.1           farver_2.1.2           
    ## INFO [2025-01-31 19:32:29]   [7] rmarkdown_2.27          GlobalOptions_0.1.2     fs_1.6.4               
    ## INFO [2025-01-31 19:32:29]  [10] zlibbioc_1.50.0         vctrs_0.6.5             memoise_2.0.1          
    ## INFO [2025-01-31 19:32:29]  [13] askpass_1.2.0           rstatix_0.7.2           htmltools_0.5.8.1      
    ## INFO [2025-01-31 19:32:29]  [16] S4Arrays_1.4.1          usethis_2.2.3           itertools_0.1-3        
    ## INFO [2025-01-31 19:32:29]  [19] missForest_1.5          lambda.r_1.2.4          curl_5.2.1             
    ## INFO [2025-01-31 19:32:29]  [22] broom_1.0.6             SparseArray_1.4.8       sass_0.4.9             
    ## INFO [2025-01-31 19:32:29]  [25] plyr_1.8.9              cachem_1.1.0            futile.options_1.0.1   
    ## INFO [2025-01-31 19:32:29]  [28] impute_1.78.0           commonmark_1.9.1        igraph_2.0.3           
    ## INFO [2025-01-31 19:32:29]  [31] lifecycle_1.0.4         iterators_1.0.14        pkgconfig_2.0.3        
    ## INFO [2025-01-31 19:32:29]  [34] Matrix_1.7-0            R6_2.5.1                fastmap_1.2.0          
    ## INFO [2025-01-31 19:32:29]  [37] GenomeInfoDbData_1.2.12 clue_0.3-65             digest_0.6.36          
    ## INFO [2025-01-31 19:32:29]  [40] pcaMethods_1.96.0       colorspace_2.1-0        RSQLite_2.3.7          
    ## INFO [2025-01-31 19:32:29]  [43] ggpubr_0.6.0            labeling_0.4.3          randomForest_4.7-1.1   
    ## INFO [2025-01-31 19:32:29]  [46] fansi_1.0.6             httr_1.4.7              abind_1.4-5            
    ## INFO [2025-01-31 19:32:29]  [49] compiler_4.4.1          rngtools_1.5.2          bit64_4.0.5            
    ## INFO [2025-01-31 19:32:29]  [52] withr_3.0.0             doParallel_1.0.17       backports_1.5.0        
    ## INFO [2025-01-31 19:32:29]  [55] DBI_1.2.3               carData_3.0-5           highr_0.11             
    ## INFO [2025-01-31 19:32:29]  [58] ggsignif_0.6.4          MASS_7.3-60.2           openssl_2.2.0          
    ## INFO [2025-01-31 19:32:29]  [61] DelayedArray_0.30.1     rjson_0.2.21            ggsci_3.2.0            
    ## INFO [2025-01-31 19:32:29]  [64] tools_4.4.1             zip_2.3.1               glue_1.7.0             
    ## INFO [2025-01-31 19:32:29]  [67] nlme_3.1-164            cluster_2.1.6           reshape2_1.4.4         
    ## INFO [2025-01-31 19:32:29]  [70] generics_0.1.3          gtable_0.3.5            tidyr_1.3.1            
    ## INFO [2025-01-31 19:32:29]  [73] xml2_1.3.6              car_3.1-2               utf8_1.2.4             
    ## INFO [2025-01-31 19:32:29]  [76] XVector_0.44.0          ggrepel_0.9.5           foreach_1.5.2          
    ## INFO [2025-01-31 19:32:29]  [79] pillar_1.9.0            markdown_1.13           stringr_1.5.1          
    ## INFO [2025-01-31 19:32:29]  [82] circlize_0.4.16         lattice_0.22-6          bit_4.0.5              
    ## INFO [2025-01-31 19:32:29]  [85] AlgDesign_1.2.1.1       tidyselect_1.2.1        juicyjuice_0.1.0       
    ## INFO [2025-01-31 19:32:29]  [88] knitr_1.47              V8_6.0.0                xfun_0.45              
    ## INFO [2025-01-31 19:32:29]  [91] credentials_2.0.1       stringi_1.8.4           UCSC.utils_1.0.0       
    ## INFO [2025-01-31 19:32:29]  [94] yaml_2.3.8              evaluate_0.24.0         codetools_0.2-20       
    ## INFO [2025-01-31 19:32:29]  [97] tibble_3.2.1            cli_3.6.3               munsell_0.5.1          
    ## INFO [2025-01-31 19:32:29] [100] Rcpp_1.0.12             gert_2.0.1              png_0.1-8              
    ## INFO [2025-01-31 19:32:29] [103] parallel_4.4.1          blob_1.2.4              doRNG_1.8.6            
    ## INFO [2025-01-31 19:32:29] [106] viridisLite_0.4.2       scales_1.3.0            openxlsx_4.2.5.2       
    ## INFO [2025-01-31 19:32:29] [109] purrr_1.0.2             crayon_1.5.3            clisymbols_1.2.0       
    ## INFO [2025-01-31 19:32:29] [112] GetoptLong_1.0.5        rlang_1.1.4             formatR_1.14
