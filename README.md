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
**iMOSCATO** requires two types of input data:

1. scRNA-seq count data, along with meta information indicating the cell type information for each cell.
2. Spatial transcriptomics count data, along with spatial location information.
   
```r
source("R/imoscato.R")
load("data/demo.Rdata")
```

### Create iMOSCATO object
The iMOSCATO object is created by the `create.iMOSCATO` function. The essential inputs are:

- sc_count: a matrix of raw scRNA-seq count data, each row represents a cell and each column represents a gene. This sc_count data serves as a reference for the cell type deconvolution for spatial transcriptomics data.
- sc_meta: a data frame of scRNA-seq meta data. The sc_meta data must contain the column indicating the cell type assignment for each cell (e.g., `cellType` column in the example sc_meta data).
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
iMOSCATO.object <- create.iMOSCATO(
  sc_count = sc_count, 
  sc_meta = sc_meta, 
  st_count = st_count, 
  loc = loc,
  cutoff_sample = 100, 
  cutoff_feature = 0.1,
  norm_method = "tss", 
  findHVG = FALSE, 
  n.HVGs=2000, 
  downsampling = FALSE, 
  size = 0.1, 
  platform = "ST")

## QC on scRNA-seq data! 
## QC on ST data! 
## Merge scRNA-seq data and ST data by common genes! 
## Joint calculation of sample-specific size factor! 
## Construct neighbor structure using ST geospatial profile!
```

### Run iMOSCATO
We run iMOSCATO using `run.iMOSCATO` function. The essential inputs are:

- iMOSCATO.object: iMOSCATO object created by `create.iMOSCATO` function.
- z: a vector of non-negative integers indicating the initial assignments of spatial domains. Starting from 0 and ending by D-1, where D is the specified number of spatial domains. It's can be obtained via K-means clustering or other clustering methods.
- iter: a number indicating the total number of iterations. Default is `5000`.
- burn: a number indicating the number of burn-in iterations. Default is `2500`.

```r
D = 2
set.seed(12345)
z = as.numeric(kmeans(iMOSCATO.object@st_count, centers = D)$cluster)-1

iMOSCATO.object = run.iMOSCATO(
  iMOSCATO.object = iMOSCATO.object, 
  z = z, 
  sc_nb = TRUE, 
  st_nb = TRUE, 
  find_domain = TRUE,
  iter = 2000,
  burn = 1000)

## iMOSCATO starts! 
10% has been done
20% has been done
30% has been done
40% has been done
50% has been done
60% has been done
70% has been done
80% has been done
90% has been done
100% has been done
## iMOSCATO finishs! Run time is 86 seconds!
```

### Deconvolute cell types
The estimated cell type proportions is stored in `iMOSCATO.object@proportion`.

```r
prop <- iMOSCATO.object@proportion
colnames(prop) <- paste0("Cell type ", 1:4)
head(prop)
              Cell type 1 Cell type 2 Cell type 3  Cell type 4
16.92x9.015     0.9699983 0.009222648  0.01496042 0.0058185899
16.945x11.075   0.1586566 0.292325599  0.53811212 0.0109056881
16.97x10.118    0.8999953 0.051818286  0.04767275 0.0005136645
16.939x12.132   0.8109649 0.017989981  0.08664395 0.0844012134
16.949x13.055   0.6669310 0.125455153  0.19360721 0.0140066167
16.942x15.088   0.9084850 0.028643868  0.06050811 0.0023630713
```
We can visualize the cell type proportion matrix through scatterpie plot via `CARD.visualize.pie` function in R package `CARD`.
```r
p1 = CARD::CARD.visualize.pie(proportion = iMOSCATO.object@proportion, 
                   spatial_location = iMOSCATO.object@loc, 
                   colors = colors) +
  theme(legend.title=element_text(size = 8),
        legend.text=element_text(size = 8),
        legend.box.spacing = unit(0, "pt"))
print(p1)
```
<img src="figure/imoscato_prop.png" alt="prop" width="325" height="300">

### Identify spatial domains
The estimated spatial domains is stored in `iMOSCATO.object@domain`.

```r
domain = iMOSCATO.object@domain
loc = iMOSCATO.object@loc
data = cbind(loc, domain)

head(data)
                    x      y domain
16.92x9.015   16.920  9.015      1
16.945x11.075 16.945 11.075      1
16.97x10.118  16.970 10.118      1
16.939x12.132 16.939 12.132      1
16.949x13.055 16.949 13.055      1
16.942x15.088 16.942 15.088      1

p2 = plot.cluster(data, x, y, size = 2, domain = as.factor(domain), colors = c("red", "steelblue3"))
print(p2)
```
<img src="figure/imoscato_domain.png" alt="domain" width="400" height="240">

### Detect discriminating genes
To obtain discriminating genes, we can check their marginal posterior probabilities of inclusion (PPI). Then, the discriminating genes are identified
if their PPI values exceed a given threshold $c$, $c$ = 0.5, which is commonly referred to as the median model.

The estimated PPI is stored in `iMOSCATO.object@PPI`.

```r
## Identified discriminating genes using c = 0.5
PPI = iMOSCATO.object@PPI
res = iMOSCATO.object@mcmc_result
feature = data.frame("gene" = colnames(st_count), "PPI" = PPI)

head(feature[feature$PPI >= 0.5, ])
     gene PPI
1   Gene1   1
4   Gene4   1
5   Gene5   1
6   Gene6   1
9   Gene9   1
19 Gene19   1

## Number of detected discriminting genes
sum(PPI >= 0.5)
[1] 10
```
