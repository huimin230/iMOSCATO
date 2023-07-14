################################################################################
## 
## Paper        : Bayesian Integrative Modeling of Single Cell and Spatial Transcriptomics Data
## Project      : iMOSCATO 
## Research goal: Bayesian cell type deconvolution with feature selection
## Data         : Mouse olfactory bulb ST data
## Modified     : 2023-07-12
## 
################################################################################

rm(list = ls())


## Load data ###################################################################
setwd("/Users/Huimin Li/PhD/Project/iMOSCATO/Code")
source("deconvolution_other_method.R")
source("functions.R")

data_name <- "mouse_olfactory_bulb"
link <- paste0("/Users/Huimin Li/PhD/Project/iMOSCATO/Data/real_data/", data_name, "/")
sc_color <- get_color("MOB", "scRNA-seq")

## Load scRNA-seq data (use normalized data instead of count data)
load(paste0(link, "processed_data/", data_name, "_scRNA_normalized.RData"))
X <- count
sc_meta <- metadata
sc_meta$cellType <- factor(sc_meta$cellType)
cellType <- levels(sc_meta$cellType)

## scRNA-seq data quality control
sc_qc <- quality_controller("scRNA-seq", count = X, metadata = sc_meta)
X <- sc_qc$count
sc_meta <- sc_qc$metadata
cell_type <- sc_meta$cellType

## Load ST data
load(paste0(link, "processed_data/", data_name, "_replicate_12.RData"))
Y <- count
loc <- loc
st_meta <- metadata
st_meta$Layer <- factor(st_meta$Layer)

## SRT data quality control
st_qc <- quality_controller("SRT", count = Y, loc = loc)
Y <- st_qc$count
loc <- st_qc$loc
Layer <- levels(st_meta$Layer)

## Common genes
sc_gene <- colnames(X)
st_gene <- colnames(Y)
common_gene <- sc_gene[sc_gene %in% st_gene]
length(common_gene)


## Running cell type deconvolution methods #####################################
methods <- c("CARD", "RCTD")

for (i in 1:2) {
  method <- methods[i]
  result <- deconvolution(method=method, X=X, Y=Y, loc=loc, sc_meta=sc_meta, st_meta=st_meta)
  
  results <- result$results
  prop_data <- result$prop_data
  domain_data <- result$domain_data
  marker_gene_data <- result$marker_gene_data
  output <- result$output
  scatterpie <- result$scatterpie
  
  
  ## Export results ##############################################################
  ## Visualization of cell type deconvolution 
  jpeg(paste0(link, "other_method/result/figure/", data_name, "_", method, "_prop.jpg"), width = 1150, height = 1000, res = 300)
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
  message(paste0("---- Method = ", method, "; ARI = ", round(output$value[1], digits = 3), "; Run time = ", round(output$value[2], digits = 0), " seconds ----"))
  
}

