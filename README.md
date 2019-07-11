
scACE
=====

Single-Cell chromatin Accessibility and gene Expression-based Clustering. A model-based approach that is specifically designed for single-cell genomic data and can jointly cluster single-cell chromatin accessibility and single-cell gene expression data.

Installation
------------

You can install the released version of scACE from Github:

``` r
library(devtools)
devtools::install_github("WWJiaxuan/scACE")
```

Main Functions
--------------

`getClusterGibbs`: Perform model-based clustering algorithm on single-cell genomic data, jointly clustering single-cell chromatin accessibility and single-cell gene expression data.

`simData`: Simulate single-cell genomic data by model-based approach, including single-cell chromatin accessibility and single-cell gene expression data for 2 clusters.

Example
-------

Please refer to the [vigenette](https://github.com/Van1yu3/SWKM/tree/master/doc) with two examples for a quick guide to scACE package.
