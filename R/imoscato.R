################################################################################
## 
## Title        : Bayesian Integrative Modeling of Single Cell and Spatial Transcriptomics Data by iMOSCATO
## Research goal: Bayesian cell type deconvolution with feature selection
## Modified     : 2023-07-13
## 
################################################################################

library(SingleCellExperiment)
library(Seurat)
library(igraph)
library(edgeR)
library(scuttle)
library(scran)
library(mcclust)
library(ggplot2)

setwd("/Users/Huimin Li/PhD/Project/iMOSCATO/Code")
source("functions.R")
Rcpp::sourceCpp("imoscato.cpp")

## Function to create an iMOSCATO object #######################################
setClass("iMOSCATO", 
         slots = list(
         sc_count = "ANY",
         sc_meta = "data.frame",
         st_count = "ANY",
         loc = "data.frame",
         v = "ANY",
         s = "ANY",
         cell_type = "ANY",
         G = "ANY",
         nei = "ANY",
         proportion = "ANY",
         domain = "ANY",
         PPM = "ANY",
         PPI = "ANY",
         mcmc_result = "list")
)


create.iMOSCATO <- function(sc_count, sc_meta, st_count, loc, cutoff_sample = 100, cutoff_feature = 0.1, 
                            norm_method = "tss", findHVG = FALSE, n.HVGs=2000, 
                            downsampling = TRUE, size = 0.1, platform = "ST"){
  
  ## Quality control (QC) ######################################################
  ## QC on scRNA-seq data 
  cat(paste0("## QC on scRNA-seq data! \n"))
  
  sc_qc <- quality_controller(data = "scRNA-seq", count = sc_count, metadata = sc_meta)
  sc_count <- sc_qc$count
  sc_meta <- sc_qc$metadata
  
  ## Check data order should be matched
  if(!identical(rownames(sc_count), rownames(sc_meta))){
    stop("The row names of sc_count and sc_meta should be should be matched each other!")
  }
  
  ## QC on ST data
  cat(paste0("## QC on ST data! \n"))
  
  st_qc <- quality_controller(data = "ST", count = st_count, loc = loc, cutoff_sample = cutoff_sample, cutoff_feature = cutoff_feature)
  st_count <- st_qc$count
  loc <- st_qc$loc
  
  ## Check data order should be matched
  if(!identical(rownames(st_count), rownames(loc))){
    stop("The row names of st_count and loc should be should be matched each other!")
  }
  
  
  ## Merge scRNA-seq data and ST data by common genes ##########################
  cat(paste0("## Merge scRNA-seq data and ST data by common genes! \n"))
  
  common_gene <- colnames(sc_count)[colnames(sc_count) %in% colnames(st_count)]
  length(common_gene)
  
  sc_count <- sc_count[, common_gene]
  st_count <- st_count[, common_gene]
  
  ## Check genes order should be matched
  check <- identical(colnames(sc_count), colnames(st_count))
  if (!check ){
    stop("The gene names of scRNA-seq data and ST data should be matched each other!")
  }
  
  
  ## Joint calculation of size factor ##########################################
  cat(paste0("## Joint calculation of sample-specific size factor! \n"))
  
  size.factor = normalize.joint(sc_count = sc_count, st_count = st_count, norm_method = norm_method)
  s = size.factor$s
  v = size.factor$v
  
  
  ## Find highly variable genes ################################################
  if (findHVG) {
    cat(paste0("## Find highly variable genes! \n"))
    
    ## Create SingleCellExperiment data using scRNA-seq data
    set.seed(1)
    rowData <- data.frame("gene" = colnames(sc_count))
    sce <- SingleCellExperiment(assays = list(counts = as(t(count), 'dgCMatrix')), rowData = rowData, colData = colData)
    sce <- logNormCounts(sce)
    dec <- modelGeneVar(sce, assay.type="logcounts")
    top <- getTopHVGs(dec, n=n.HVGs)
    sc_count <- sc_count[, top]
    st_count <- st_count[, top]
    
  }
  
  ## Cell downsampling of scRNA-seq data #######################################
  if (downsampling) {
    cat(paste0("## Cell downsampling of scRNA-seq data! \n"))
    
    ## Create Seurat object and run dimension reduction
    seurat_object = create_seurat(sc_count, sc_meta)

    ## Cell downsampling
    selected_cell = down_sampling(seurat_object, size = size)

    sc_count <- sc_count[selected_cell, ]
    sc_meta <- sc_meta[selected_cell, ]
    v <- v[selected_cell]
    
  }
  
  
  ## Cell type information #####################################################
  sc_meta$cellType <- factor(sc_meta$cellType)
  cellType <- levels(sc_meta$cellType)
  K <- length(cellType)
  cell_type <- as.numeric(factor(sc_meta$cellType, levels = cellType, labels = 1:K)) - 1
  
  
  ## Construct neighbor structure ##############################################
  cat(paste0("## Construct neighbor structure using ST geospatial profile! \n"))
  
  if (platform == "ST"){
    tt <- get.neighbor(loc, 4)
    dist <- tt$dist
    G <- tt$G
    nei <- tt$nei
  }else if(platform == "Visium"){
    tt <- get.neighbor(loc, 6)
    dist <- tt$dist
    G <- tt$G
    nei <- tt$nei
  }else{
    ## Voronoi tessellation
    temp <- data.frame(id = 1:n, x = loc$y, y = loc$x)
    tt = voronoi_adjacency(data = temp, id~x+y, scale=1, PLOT=FALSE)
    dist <- tt$dist
    G <- tt$G
    nei <- tt$nei
    
  }
  
  pp <- which(colSums(nei) == 0)
  if (length(pp) > 0){
    P <- P[, -pp]
  }
  

  ## Create and export the object ##############################################
  object <- new(Class = "iMOSCATO",
                sc_count = sc_count,
                sc_meta = sc_meta,
                st_count = st_count,
                loc = loc,
                v = v,
                s = s,
                cell_type = cell_type,
                G = G,
                nei = nei)
  
  return(object)

}


## Function to run iMOSCATO ####################################################
run.iMOSCATO <- function(iMOSCATO.object, z, sc_nb = TRUE, st_nb = TRUE, find_domain = TRUE, iter = 2000, burn = 1000){
  
  sc_count = iMOSCATO.object@sc_count
  sc_meta = iMOSCATO.object@sc_meta
  st_count = iMOSCATO.object@st_count
  loc = iMOSCATO.object@loc
  v = iMOSCATO.object@v
  s = iMOSCATO.object@s
  cell_type = iMOSCATO.object@cell_type
  G = iMOSCATO.object@G
  nei = iMOSCATO.object@nei
  
  sc_meta$cellType <- factor(sc_meta$cellType)
  cellType <- levels(sc_meta$cellType)
  
  cat(paste0("## iMOSCATO starts! \n"))
  start_time = proc.time()
  
  res = imoscato(X=as.matrix(sc_count), Y=as.matrix(st_count), cell_type=cell_type, 
                 z=z, v=v, s=s, P=nei, sc_nb=sc_nb, st_nb=st_nb, find_domain=find_domain, iter=iter, burn=burn)
  
  end_time = proc.time()
  time = as.numeric((end_time - start_time)[3], "secs")
  cat(paste0("## iMOSCATO finishs! Run time is ", round(time, digits = 0), " seconds! \n"))
  
  object = iMOSCATO.object
  res$time = time
  object@mcmc_result = res
  
  
  ## Estimated cell type proportions
  Pi_store = lapply(1:iter, function(x) res$Pi_store[ , , x])
  proportion = Reduce("+", Pi_store[(burn + 1):iter]) / (iter - burn)
  
  rownames(proportion) = rownames(loc)
  colnames(proportion) = cellType
  object@proportion = proportion
  
  ## Estimated spatial domain
  D = length(unique(z))
  tt = get.ppm(res$z_store, burn, iter, D)
  domain = tt$z_ppm
  object@domain = domain
  PPM = tt$ppm
  rownames(PPM) = rownames(loc)
  colnames(PPM) = rownames(loc)
  object@PPM = PPM
  
  ## Identify discriminating genes
  object@PPI = res$gamma_ppi
  
  return(object)
  
}
