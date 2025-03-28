################################################################################
## 
## Title        : Robust Bayesian Integrative Modeling of Single Cell and Spatially Resolved Transcriptomics Data
## Project      : iMOSCATO
## Research goal: Bayesian cell-type deconvolution with spatial domain detection
## Authors      : Huimin Li
## Contact      : huimin.li@utdallas.edu
## Modified     : 2025-04-01
## 
################################################################################


library(Rcpp)
library(SingleCellExperiment)
library(Seurat)
library(igraph)
library(edgeR)
library(scuttle)
library(scran)
library(mcclust)
library(ggplot2)
library(gridExtra)
library(ggpubr)

source("R/functions.R")
source("R/deconvolution_method.R")
sourceCpp("R/imoscato.cpp") 
which <- base::which


## Function to create an iMOSCATO object #######################################
setClass("iMOSCATO", 
         slots = list(
           sc_count = "ANY",
           sc_meta = "data.frame",
           B = "ANY",
           st_count = "ANY",
           loc = "data.frame",
           v = "ANY",
           s = "ANY",
           G = "ANY",
           nei = "ANY")
)


create.iMOSCATO <- function(sc_count, sc_meta, st_count, loc, cutoff_sample = 100, cutoff_feature = 0.1, 
                            norm_method = "tss", platform = "ST"){
  
  
  ## Quality control (QC) ######################################################
  # QC on scRNA-seq data
  cat(paste0("## Quality control on scRNA-seq data! \n"))

  ## Check data order should be matched
  if(!identical(rownames(sc_count), rownames(sc_meta))){
    stop("The row names of count and metadata should be should be matched each other!")
  }
  
  sc_qc <- sc.qc(count = sc_count, metadata = sc_meta)
  sc_count <- sc_qc$count
  sc_meta <- sc_qc$metadata

  # print(paste0("Number of cells is ", nrow(sc_count), ", number of genes is ", ncol(sc_count), " after quality control"))
  
  ## QC on SRT data
  cat(paste0("## Quality control on spatially resolved transcriptomics data! \n"))
  
  ## Check data order should be matched
  if(!identical(rownames(st_count), rownames(loc))){
    stop("The row names of count and loc should be matched each other!")
  }
  
  st_qc <- st.qc(count = st_count, loc = loc, cutoff_sample = cutoff_sample, cutoff_feature = cutoff_feature)
  st_count <- st_qc$count
  loc <- st_qc$loc
  
  # print(paste0("Number of spots is ", nrow(st_count), ", number of genes is ", ncol(st_count), " after quality control"))
  
  
  ## Calculation of size factor ################################################
  common_gene <- colnames(st_count)[colnames(st_count) %in% colnames(sc_count)]
  # print(paste0("Number of common genes is ", length(common_gene), " after quality control"))
  
  cat(paste0("## Joint calculation of sample-specific size factor! \n"))
  
  size.factor = normalize.joint(sc_count = sc_count[, common_gene], st_count = st_count[, common_gene], norm_method = norm_method)
  s = size.factor$s
  v = size.factor$v
  
  
  ## Calculate reference basis matrix B from scRNA-seq data ####################
  cat(paste0("## Create reference basis matrix from scRNA-seq data! \n"))

  result <- quiet(CARD_run(sc_count=sc_count, sc_meta=sc_meta, st_count=st_count, loc=loc))
  B <- result$reference_matrix
  sc_count <- sc_count[, colnames(B)]
  st_count <- st_count[, colnames(B)]

  ## Check genes order should be matched
  if (!identical(colnames(sc_count), colnames(st_count))){
    stop("The gene names of scRNA-seq data and SRT data should be matched each other!")
  }
  
  
  ## Construct neighbor structure ##############################################
  cat(paste0("## Construct neighbor structure using spatial transcriptomics geospatial profile! \n"))
  
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
    nei <- nei[, -pp]
  }

  rownames(nei)=rownames(loc)  

  
  ## Create and export the object ##############################################
  object <- new(Class = "iMOSCATO",
                sc_count = sc_count,
                sc_meta = sc_meta,
                B = B,
                st_count = st_count,
                loc = loc,
                v = v,
                s = s,
                G = G,
                nei = nei)
  
  return(object)

}


## Function to run iMOSCATO ####################################################
run.iMOSCATO <- function(iMOSCATO.object, n.domain, iter = 2000, burn = 1000){
  
  B = iMOSCATO.object@B
  st_count = iMOSCATO.object@st_count
  loc = iMOSCATO.object@loc
  s = iMOSCATO.object@s
  G = iMOSCATO.object@G
  nei = iMOSCATO.object@nei

  cat(paste0("## iMOSCATO starts! \n"))
  start_time = proc.time()
  
  cellType <- rownames(B)
  
  z <- suppressWarnings(cluster_start(count=iMOSCATO.object@st_count, N_domain = n.domain))
  
  res = imoscato(Y=as.matrix(st_count), 
                 B=as.matrix(B),
                 z=z, 
                 s=s, 
                 P=nei, 
                 iter=iter, 
                 burn=burn)
  
  end_time = proc.time()
  time = as.numeric((end_time - start_time)[3], "secs")
  cat(paste0("## iMOSCATO finishs! Run time is ", round(time, digits = 0), " seconds! \n"))
  
  result <- list()
  mcmc_result = res
  object = iMOSCATO.object
  
  result$time = time

  ## Estimated cell type prop_results
  prop_result = res$Pi_map
  rownames(prop_result) = rownames(loc)
  colnames(prop_result) = cellType
  result$prop_result = as.data.frame(prop_result)
  
  ## Estimated spatial domain
  D = length(unique(z))
  tt = get.ppm(res$z_store, burn, iter, D)
  domain = tt$z_ppm
  PPM = tt$ppm
  rownames(PPM) = rownames(loc)
  colnames(PPM) = rownames(loc)
  result$PPM = PPM
  
  dominant_type <- unlist(sapply(1:dim(prop_result)[1], function(x) cellType[which.max(prop_result[x, ])]))
  
  domain_result <- data.frame("x" = loc$x, "y" = loc$y, "domain" = domain, "domain_map" = res$z_map+1, "dominant_type"= dominant_type)
  rownames(domain_result) <- rownames(loc)
  result$domain_result = domain_result

  
  ## Other parameters
  Omega_hat = res$Omega_map
  Delta_hat = ifelse(res$Delta_ppi>=0.5, 1, 0)
  A_hat = res$A_map
  phi_hat = res$phi_map
  Eta_hat = ifelse(res$H_ppi>=0.5, 1, 0)

  rownames(Omega_hat) = rownames(Delta_hat) = rownames(Eta_hat) = rownames(loc)
  colnames(Omega_hat) = colnames(Delta_hat) = colnames(A_hat) = cellType
  colnames(Eta_hat) = names(phi_hat) = colnames(st_count)

  result$Omega_hat = Omega_hat
  result$Delta_hat = Delta_hat
  result$A_hat = A_hat
  result$phi_hat = phi_hat
  result$Eta_hat = Eta_hat
  
  z_max <- which.max(table(domain_result$domain_map))-1
  pi_max <- which.max(prop_result[1,])
  
  ## Check the convergence
  trace_data <- data.frame("Iteration" = 1:iter, 
                           "pi" = res$Pi_store[1,pi_max,],
                           "delta" = res$Delta_sum,
                           "z_max" = sapply(1:iter, function(x)sum(res$z_store[x,]==z_max)))
  
  trace_fig <- c()
  
  for (i in 1:ncol(trace_data)) {
    y = trace_data[1, i]
    x = trace_data[, i]
    trace_data[, i][x > quantile(x, probs = 0.999)] <- quantile(x, probs = 0.999)
    trace_data[, i][x < quantile(x, probs = 0.001)] <- quantile(x, probs = 0.001)
    trace_data[1, i] = y
  }
  
  ## Trace plot of one prop_result
  trace_fig[[1]] <- plot.trace(x = trace_data$Iteration, y = trace_data$pi, 
                               ylab = "1st Spot largest cell type proportion",
                               ylim = c(min(trace_data$pi), 1))
  
  ## Trace plot of number of zero prop_results
  trace_fig[[2]] <- plot.trace(x = trace_data$Iteration, y = trace_data$delta, 
                               ylab = "Number of zero proportion",
                               ylim = c(min(trace_data$delta), 10 + max(trace_data$delta)))
  
  ## Trace plot of number of spots in spatial domain with most spots
  trace_fig[[3]] <- plot.trace(x = trace_data$Iteration, y = trace_data$z_max, 
                               ylab = "Number of spots in largest domain",
                               ylim = c(min(trace_data$z_max), 5 + max(trace_data$z_max)))
  
  result$trace_result <- trace_data

  return(list("object"=iMOSCATO.object,"mcmc_result"=mcmc_result,"result"=result,"trace_fig"=trace_fig))
  
}


