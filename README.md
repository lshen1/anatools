# anatools

[![Build Status](https://travis-ci.org/lshen1/anatools.svg?branch=master)](https://travis-ci.org/lshen1/anatools)

The **anatools** R package is a collection of analytic tools for analyzing high-throughput data,
including microarray and next-generation-sequencing to quickly produce results, tables and charts.
It can be used to:

+ utility_function:
  * To show, get, rename an R object.
  * To aggregate multiple probes of microarry data based on the mean or maximun variance of gene probesets.
+ modelTesting:
  * Perform Two-Sample t-tests, one-way ANOVA with contrasts using `limma` package in R.
  * Perform one-way ANOVA with TukeyHSD using `aov` package in R.
  * Perform Negative binomial generalized linear models using `DEseq2` package in R.
  * Perform Correlation tests (pearson/spearman) using `Hmisc` package in R.
  * Perform univariate cox model for expression data.
  * Perform Fisher's exact test and Pearson's Chi-squared Test using `stats` package in R.
  * Performs Wilcoxon tests using `stats` package in R.
+ autoplot:
  * Generate Kaplan-Meier Plot and table using `survMisc` package in R.
+ BI:
  * Generate density-plot with Bimodality Index (BI) for microarray and separate the groups based on BI.
+ color:
  * To diverge colour palette function with set midpoint.
  * Add transparent/alpha to a set of colours.
+ ggplot_wrapper:
  * Generate Heatmap using `NMF` package in R.
  * Generate a plot of beta-uniform mixture model to a set of p-values (BUM plot) with significant counts table based on a set of FDRs.
  * Generate a plot of table to show the features whose coefficient greater than cutoffs.
  * Generate a plot of table to show the results of Fisher's Exact Test for Count Data.
  * Generate density-plot for QC accessment on microarray data.
  * Generate box-plot for QC accessment on microarray data.
  * Generate principal components analysis (PCA) Plot on microarray data.
  * Generate correlation-plot between X and Y.
  * Generate concordance-plot between X and Y.
  * Generate bar-plot for the frequency of elements in data.frame.
  * Make a Venn Diagram (VD) with count table that each component corresponding to a separate circle in VD.
  * Generate a Quantile-Quantile Plot (Q-Q plot).

## Installation

The [**devtools** package](http://cran.r-project.org/web/packages/devtools/index.html) is used to install R packages hosted on Github. To install **anatools**, type the following commands in the R console:
```r
    library(devtools)
    install_github("lshen1/anatools")
```

## Usage

```r
# load the package
library("anatools")
```








