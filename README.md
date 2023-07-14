# iMOSCATO

Bayesian <ins>i</ins>ntergrative <ins>M</ins><ins>O</ins>deling of <ins>S</ins>ingle <ins>C</ins>ell and Sp<ins>A</ins>tial <ins>T</ins>ranscript<ins>O</ins>mics Data

## Introduction

**iMOSCATO** is a fully Bayesian model that integrates spatial transcriptomics data and single-cell RNA sequencing (scRNA-seq) data to decompose cell-type mixtures of regularly distributed spots and identify the underlying spatial domains simultaneously. It incorporates the lattice structure by employing a Markov random field prior to improve the accuracy of cell-type deconvolution and spatial domain detection. Moreover, iMOSCATO employ a feature selection mechanism to improve the algorithm efficiency, while allowing the discovery of discriminating genes whose expression levels are significantly different among cell types.

![iMOSCATO](figure/imoscato_workflow.png)

**iMOSCATO** was developed and tested under `R 4.2.3`. The following R packages are required to run the model

- Rcpp
- RcppArmadillo
- RcppDist
- SingleCellExperiment
- Seurat
- igraph
- edgeR
- scuttle
- scran
- mcclust

## Run iMOSCATO on demo data

The following section will guide to run a exemplary data using **iMOSCATO**.

### Load iMOSCATO and demo data
```r
source("R/imoscato.R")
load("data/demo.Rdata")
```

### Create an iMOSCATO object
The iMOSCATO object is created by the function `create.iMOSCATO`. The essential inputs are:

- sc_count: a matrix of raw scRNA-seq count data, each row represents a cell and each column represents a gene. This sc_count data serves as a reference for the cell type deconvolution for spatial transcriptomics data.
- sc_meta: a data frame of scRNA-seq metadata data. The sc_meta data must contain the column indicating the cell type assignment for each cell (e.g., `cellType` column in the example sc_meta data).
- st_count: a matrix of raw spatial transcriptomics count data, each row represents a spot and each column represents a gene. This is the spatial transcriptomics data that we are interested to deconvolute.
- loc: a data frame with two columns representing the $x$ and $y$ coordinates of the spot.
- cutoff_sample: a number indicating that spots are kept with at least this number of total counts across all genes. Default is 100.
- cutoff_feature: a number indicating that genes are kept with at least this percent of spots with non-zero counts. Default is 0.1.
- norm_method: a character string specifying the method to calculate the sample-specific size factor, must be one of `tss`, `q75`, `rle`, or `tmm`. Default is `tss`.
- findHVG: a logical variable indicating whether to find the highly variable genes. Default is `FALSE`.
- n.HVGs: a number indicating the number of highly variable genes to be selected. Default is `2000`.
- downsampling: a logical variable indicating whether to perform cell downsampling of scRNA-seq data. Default is `TRUE`.
- size: a way to select the representative cells in downsampling, must be a percent (select same percent for each cell type), a positive integer (select certain number of cells, then each selected cell type will have the similar percent as the orginal data), a percent vector (select desired percent for each cell type), or a positive integer vector (select desired number of cells for each cell type). Default is `0.1`.
- platform: a character string specifying the ST technology in order to construct neighbor structure, must be one of `ST`, `Visium`, or `other` (for any technologies other than `ST` and `10x Visium`). Default is `ST`.


```r
result <- dataPreprocess(
  count = count, 
  loc = loc, 
  cutoff_sample = 100, 
  cutoff_feature = 0.1, 
  cutoff_max = 0, 
  size.factor = "tss", 
  platform = "ST",
  findHVG = FALSE, 
  n.HVGs=2000)

count <- result$count
loc <- result$loc
s <- result$s
P <- result$P
```

### Run the model
We run the model using function `bayes_cafe`, where `K` is the specified number of clusters.

```r
res <- bayes_cafe(
  count = count, 
  loc = loc, 
  K = 2, 
  s = s, 
  P = P)
```

### Decompose cell type

### Identify the discriminating genes
The main purpose of **BayesCafe** is to identify discriminating genes and cluster spatial locations.
To obtain discriminating genes, we can check their marginal posterior probabilities
of inclusion (PPI). Then, the discriminating genes are identified
if their PPI values exceed a given threshold $c$, $c$ = 0.5, which is commonly referred to as the median model.

Alternatively, we can determine the threshold that controls for multiplicity, which ensures that the expected Bayesian false discovery rate (BFDR) is less than a
specified value (e.g. 0.05). 

```r
## Identified discriminating genes using c = 0.5
head(res$gamma[res$gamma$PPI >= 0.5, ])
      gene PPI
7   gene 7   1
17 gene 17   1
20 gene 20   1
36 gene 36   1
40 gene 40   1
46 gene 46   1

sum(res$gamma$PPI >= 0.5)
[1] 15

## Identified discriminating genes to control BFDR < 0.05
(threshod <- bfdr(PPI = res$gamma$PPI, alpha = 0.05))
[1] 0.9

head(res$gamma[res$gamma$PPI > 1 - threshod, ])
      gene PPI
7   gene 7   1
17 gene 17   1
20 gene 20   1
36 gene 36   1
40 gene 40   1
46 gene 46   1

sum(res$gamma$PPI > 1 - threshod)
[1] 15
```


### Visualize the clustering results
```r
head(res$cluster)
                   x      y     cluster
16.92 x 9.015   16.920  9.015       1
16.945 x 11.075 16.945 11.075       1
16.97 x 10.118  16.970 10.118       1
16.939 x 12.132 16.939 12.132       1
16.949 x 13.055 16.949 13.055       1
16.942 x 15.088 16.942 15.088       1

plot.cluster(res$cluster, x, y, group = as.factor(cluster), colors = c("red", "steelblue3"))
```
<img src="cluster.png" alt="cluster" width="500" height="300">
