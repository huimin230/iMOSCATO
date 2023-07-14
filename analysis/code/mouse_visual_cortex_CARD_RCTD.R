################################################################################
## 
## Paper        : Bayesian Integrative Modeling of Single Cell and Spatial Transcriptomics Data
## Project      : iMOSCATO 
## Research goal: Bayesian cell type deconvolution with feature selection
## Data         : Mouse visual cortex STARmap data
## Modified     : 2023-07-12
## 
################################################################################

rm(list = ls())

setwd("/Users/Huimin Li/PhD/Project/iMOSCATO/Code")
source("deconvolution_other_method.R")
source("functions.R")

data_name <- "mouse_visual_cortex"
link <- paste0("/Users/Huimin Li/PhD/Project/iMOSCATO/Data/real_data/", data_name, "/")
sc_color <- get_color("STARmap", "scRNA-seq")
library(SingleCellExperiment)


## Load data ###################################################################
## Load scRNA-seq data
load(paste0(link, "processed_data/", data_name, "_scRNA.RData"))
X <- count
sc_meta <- metadata

## scRNA-seq data quality control
sc_qc <- quality_controller("scRNA-seq", count = X, metadata = sc_meta)
X <- sc_qc$count
sc_gene <- colnames(X)
sc_meta <- sc_qc$metadata
sc_meta$cellType <- factor(sc_meta$cellType)
cellType <- levels(sc_meta$cellType)

## Log normalization
counts <- t(X)
rowData <- data.frame("gene" = rownames(counts))
colData <- sc_meta

sce <- SingleCellExperiment(assays = list(counts = as(counts, 'dgCMatrix')), rowData = rowData, colData = colData)
sce <- scuttle::logNormCounts(sce)
sc_norm <- as.data.frame(t(as.matrix(sce@assays@data@listData$logcounts)))


## Load STARmap spot level data
load(paste0(link, "processed_data/", data_name, "_STARmap_spot_level_grid_750.RData"))
Y <- count
loc <- loc
st_meta <- metadata
st_meta$Layer <- factor(st_meta$Layer)
true_prop <- prop

## SRT data quality control
st_qc <- quality_controller("SRT", count = Y, loc = loc)
Y <- st_qc$count
loc <- st_qc$loc

Layer <- levels(st_meta$Layer)
st_gene <- colnames(Y)

## Common gene
common_gene <- sc_gene[sc_gene %in% st_gene]
length(common_gene)


## Running cell type deconvolution methods #####################################
methods <- c("CARD", "RCTD")

for (i in 1:2) {
  method <- methods[i]
  
  if (method == "CARD"){
    X = sc_norm
  } else{
    X = X
  }
  
  result <- deconvolution(method=method, X=X, Y=Y, loc=loc, sc_meta=sc_meta, st_meta=st_meta, true_prop=true_prop)
  
  results <- result$results
  prop_data <- result$prop_data
  marker_gene_data <- result$marker_gene_data
  output <- result$output
  
  spatial_location = prop_data[, c("y", "x")]
  names(spatial_location) <- c("x", "y")
  scatterpie <- CARD.visualize.pie(proportion = prop_data[, cellType], spatial_location = spatial_location, colors = sc_color)

  
  ## Export results ##############################################################
  ## Visualization of cell type deconvolution 
  jpeg(paste0(link, "other_method/result/figure/", data_name, "_", method, "_prop.jpg"), width = 1300, height = 800, res = 300)
  gridExtra::grid.arrange(scatterpie, nrow = 1)
  dev.off()
  
  ## Export results
  save(results, file = paste0(link, "other_method/output/", data_name, "_", method, "_result.RData"))
  
  ## Export estimated cell type proportion 
  write.csv(prop_data, paste0(link, "other_method/result/", data_name, "_", method, "_prop.csv"), row.names = FALSE)
  
  ## Export discriminating genes
  write.csv(marker_gene_data, paste0(link, "other_method/result/", data_name, "_", method, "_marker_gene.csv"), row.names = FALSE)
  
  ## Export output results
  write.csv(output, paste0(link, "other_method/result/",  data_name, "_", method, "_metrics.csv"), row.names = FALSE)
  
  ## Print results
  message(paste0("---- Method = ", method, "; RMSE = ", round(output$value[1], digits = 3), "; PCC = ",  round(output$value[2], digits = 3), "; ARI = ", round(output$value[3], digits = 3), "; Run time = ", round(output$value[4], digits = 0), " seconds ----"))
  
}

