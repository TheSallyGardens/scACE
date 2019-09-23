
scACE
=====

Single-Cell chromatin Accessibility and gene Expression-based Clustering. A model-based approach that is specifically designed for single-cell genomic data and can jointly cluster single-cell chromatin accessibility and single-cell gene expression data.

![](https://github.com/WWJiaxuan/scACE/blob/master/scACE_pic.jpg)

Installation
------------

You can install the released version of scACE from Github:

``` r
library(devtools)
devtools::install_github("cuhklinlab/scACE")
```

Main Functions
--------------

`getClusterGibbs`: Perform model-based clustering algorithm on single-cell genomic data using Markov Chain Monte Carlo (MCMC), jointly clustering single-cell chromatin accessibility and single-cell gene expression data.

`update_all2`: Perform model-based clustering algorithm on single-cell genomic data using expectation–maximization (EM) algorithm, jointly clustering sc-ATAC and sc-RNA data.

`simData`: Simulate single-cell genomic data by model-based approach, including single-cell chromatin accessibility and single-cell gene expression data for 2 clusters.

Example
-------

Please refer to the [vigenette](https://github.com/WWJiaxuan/scACE/tree/master/vignette) with several examples for a quick guide to scACE package.
