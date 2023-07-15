################################################################################
## 
## Title        : Bayesian Integrative Modeling of Single Cell and Spatial Transcriptomics Data by iMOSCATO
## Research goal: Bayesian cell type deconvolution with feature selection
## Modified     : 2023-07-13
## 
################################################################################


## Data quality control ########################################################
quality_controller <- function(data, count, loc = NULL, metadata = NULL, cutoff_sample = 100, cutoff_feature = 0.1) {
  if (data == "scRNA-seq"){
    ## QC on scRNA-seq data
    m <- dim(count)[1]
    p <- dim(count)[2]
    
    ## Sample-wise quality control for scRNA-seq data
    index <- which(rowSums(count > 0) > 0)
    count <- count[index, ]
    metadata <- metadata[index, ]
    
    ## Feature-wise quality control for scRNA-seq data
    index <- which(colSums(count > 0) > 0)
    count <- count[, index]
    
    return (list("count" = count, "metadata" = metadata))
    
  }else{
    ## QC on ST data
    n <- dim(count)[1]
    p <- dim(count)[2]
    
    ## Sample-wise quality control for ST data
    index <- which(rowSums(count) >= cutoff_sample)
    count <- count[index, ]
    loc <- loc[index, ]
    
    ## Feature-wise quality control for ST data
    index <- which(colSums(count > 0) >= cutoff_feature*n)
    count <- count[, index]
    
    return (list("count" = count, "loc" = loc))
  }
}


## Create Seurat objects and run dimension reduction ###########################
create_seurat <- function(count, metadata){
  ## Construct Seurat object 
  set.seed(0)
  K = length(unique(metadata$cellType))
  gene_name = data.frame("gene" = colnames(count))
  
  seurat_object = CreateSeuratObject(counts = t(count), project = "seurat_object", assay = "RNA", 
                                     meta.data = metadata, meta.features = gene_name, 
                                     min.cells = 0, min.features = 0)
  rownames(gene_name) = gene_name$gene
  seurat_object[["RNA"]]@meta.features = gene_name[gene_name$gene %in% rownames(seurat_object), -1]
  
  ### filtering
  seurat_object <- subset(seurat_object, subset = nFeature_RNA > 100)
  
  ### normalizing data
  seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
  
  ### find variable genes
  seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 500)
  
  ### scale the data before dimension reduction
  all.genes <- rownames(seurat_object)
  seurat_object <- ScaleData(seurat_object, features = all.genes)
  
  ### run PCA with the variable genes
  seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))
  
  #### run other dimension reduction algorithm
  seurat_object <- RunTSNE(seurat_object, features = VariableFeatures(object = seurat_object))
  seurat_object <- RunUMAP(seurat_object, features = VariableFeatures(object = seurat_object), n.components = 15)
  
  return(seurat_object)
}


## Cell downsampling ###########################################################
## size could be the number of selected cells or the proportion
down_sampling <- function(seurat_object, size = 0.1){
  seurat_object = FindNeighbors(seurat_object, features = VariableFeatures(object = seurat_object), k.param = 10)
  
  cell_graph = graph_from_adjacency_matrix(seurat_object@graphs$RNA_snn, mode = "undirected", weighted = "weight")
  
  ### construct subgraphs for each cell type
  cell_prop = table(seurat_object@meta.data$cellType) / dim(seurat_object)[2]
  
  if (length(size) == 1){
    if(size < 1 & size > 0){
      cells = ceiling(cell_prop * size *dim(seurat_object)[2])
    } else if (size >= 1){
      cells = ceiling(cell_prop * size)
    } else{
      stop("size should be a positive integer!")
    }
    
  } else if (length(size) == K) {
    if(max(size) <= 1 & min(size) > 0){
      cells = ceiling(cell_prop * size *dim(seurat_object)[2])
    } else if (min(size) > 1){
      cells = size
      names(cells) = names(cell_prop)
    } else{
      stop("size should be a positive integer vector of the length same as the number of cell types!")
    }
    
  } else{
    stop("size should be a positive integer or a positive integer vector of the length same as the number of cell types!")
  }
  
  ### sub_graph of each cell type
  selected_cell = c()
  for(k in 1:length(cells)){
    type_graph = subgraph(cell_graph, vids = which(seurat_object@meta.data$cellType == names(cells)[k]))
    result = page.rank(type_graph, algo = "prpack", directed = F)
    selected_cell_k = names(result$vector[order(result$vector, decreasing = TRUE)][1:cells[k]])
    selected_cell = c(selected_cell, selected_cell_k)
  }
  
  return(selected_cell)
}


## st data normalization #######################################################
normalize.st <- function(count, norm_method = 'tss'){
  library(edgeR)
  gene_num <- ncol(count)
  sample_num <- nrow(count)
  
  count_rowsum <- rowSums(count)
  
  if(norm_method == "tss")
  {
    ## TSS(Total Sum Scaling)
    ### scale-factors
    raw_s_factors <- rowSums(count)
    scale_coeff <- exp((-1/nrow(count)) * sum(log(raw_s_factors)))
    scaled_s_factors <- scale_coeff * raw_s_factors
    
    ### normalized count matrix
    db.norm <- sweep(count, 1, rowSums(count), FUN = "/")
    count_nor <- db.norm
  }
  else if(norm_method == "q75")
  {
    ## Q75(Upper Quantile normalization)
    ### scale factors (remove those with 0's)
    raw_s_factors <- apply(count, 1,function(x){quantile(x[x>0],0.75)} )
    scale_coeff <- exp((-1/nrow(count)) * sum(log(raw_s_factors)))
    scaled_s_factors <- scale_coeff * raw_s_factors
    
    ### normalized count matrix
    db.norm <- sweep(count, 1, raw_s_factors, FUN = "/")
    count_nor <- db.norm
  }
  else if(norm_method == "rle")
  {
    ## RLE(Relative Log Expression normalization)
    ### scale_factors
    ### function for calculating the geometric mean
    geo_mean <- function(x){
      exp(sum(log(x[x>0]))/length(x))
    }
    ### function for calculating non-zero median
    non_zero_median <- function(x){
      median(x[x>0])
    }
    ref_sample <- apply(count, 2, geo_mean)
    norm_rle_1 <- sweep(count, 2, ref_sample, FUN = "/")
    raw_s_factors <- apply(as.matrix(norm_rle_1), 1, non_zero_median)
    scale_coeff <- exp((-1/nrow(count)) * sum(log(raw_s_factors)))
    scaled_s_factors <- scale_coeff * raw_s_factors
    
    ### normalized count matrix
    db.norm <- sweep(count, 1, raw_s_factors, FUN = "/")
    count_nor <- db.norm
  }
  else if(norm_method == "tmm")
  {
    ## TMM(Trimmed Mean Method)
    ### scale_factors
    count_t <- t(count)
    raw_s_factors <- calcNormFactors(count_t, method = "TMM")
    scale_coeff <- exp((-1/nrow(count)) * sum(log(raw_s_factors)))
    scaled_s_factors <- scale_coeff * raw_s_factors
    
    # normalized count matrix
    db.norm <- sweep(count, 1, raw_s_factors, FUN = "/")
    count_nor <- db.norm
  }
  
  else if(norm_method == "n-vst") {
    ## Naive Transformation (VST)
    varx = apply(count, 2, var)
    meanx = apply(count, 2, mean)
    phi = coef(nls(varx ~ meanx + phi * meanx^2, start = list(phi = 1)))
    y <- sqrt( phi * count)
    naive_norm <- log(y + sqrt(1 + y^2))
    total_n_counts <- apply(count, 1, sum)
    log_total_n_counts <- log(total_n_counts)
    db.norm <- apply(naive_norm, 2, function(x){resid(lm(x ~ log_total_n_counts))} )
    ## All the above normalized counts were negative so reversed their signs
    count_nor <- db.norm
  }
  else if(norm_method == "a-vst")
  {
    ## Anscombe's Transformation(VST)
    varx = apply(count, 2, var)
    meanx = apply(count, 2, mean)
    phi = coef(nls(varx ~ meanx + phi * meanx^2, start = list(phi = 1)))
    ### Took the absolute value of counts as negative under the square root not on Real subspace
    y <- sqrt( abs( (count + (3/8))/((1/phi)-(3/4)) ) )
    anscombe_norm <- log(y + sqrt(1 + y^2))
    total_a_counts <- apply(count, 1, sum)
    log_total_a_counts <- log(total_a_counts)
    db.norm <- apply(anscombe_norm, 2, function(x){resid(lm(x ~ log_total_a_counts))} )
    ## All the above normalized counts were negative so reversed their signs
    count_nor <- db.norm
  }
  else if(norm_method == "log-vst")
  {
    ## Log Transformation
    varx = apply(count, 2, var)
    meanx = apply(count, 2, mean)
    phi = coef(nls(varx ~ meanx + phi * meanx^2, start = list(phi = 1)))
    log_norm <- log2(count + (1/(2*phi)))
    total_l_counts <- apply(count, 1, sum)
    log_total_l_counts <- log(total_l_counts)
    db.norm <- apply(log_norm, 2, function(x){resid(lm(x ~ log_total_l_counts))} )
    ## All the above normalized counts were negative so reversed their signs
    count_nor <- db.norm
  }
  else if(norm_method == 'none'){
    scaled_s_factors <- rep(1, sample_num)
    count_nor <- count
  }
  else if(norm_method == "none"){
    scaled_s_factors = rep(1, sample_num)
    
  }
  else
  {
    stop("Please choose a valid normalization method")
  }
  colnames(count_nor) <- colnames(count)
  return(scaled_s_factors)
}


## Data normalization jointly ##################################################
normalize.joint <- function(sc_count, st_count, norm_method = 'tss'){
  n <- nrow(st_count)
  m <- nrow(sc_count)
  
  st_count_rowsum <- rowSums(st_count)
  sc_count_rowsum <- rowSums(sc_count)
  
  if(norm_method == "tss")
  {
    ## TSS(Total Sum Scaling)
    ### scale-factors
    raw_s_factors <- rowSums(st_count)
    raw_v_factors <- rowSums(sc_count)
    
    scale_coeff <- exp((-1/(n + m)) * (sum(log(raw_s_factors)) + sum(log(raw_v_factors))))
    scaled_s_factors <- scale_coeff * raw_s_factors
    scaled_v_factors <- scale_coeff * raw_v_factors
  }
  else if(norm_method == "q75")
  {
    ## Q75(Upper Quantile normalization)
    ### scale factors (remove those with 0's)
    raw_s_factors <- apply(st_count, 1, function(x){quantile(x[x>0],0.75)} )
    raw_v_factors <- apply(sc_count, 1, function(x){quantile(x[x>0],0.75)} )
    
    scale_coeff <- exp((-1/(n + m)) * (sum(log(raw_s_factors)) + sum(log(raw_v_factors))))
    scaled_s_factors <- scale_coeff * raw_s_factors
    scaled_v_factors <- scale_coeff * raw_v_factors
  }
  else if(norm_method == "rle")
  {
    ## RLE(Relative Log Expression normalization)
    ### scale_factors
    ### function for calculating the geometric mean
    geo_mean <- function(x){
      exp(sum(log(x[x>0]))/length(x))
    }
    ### function for calculating non-zero median
    non_zero_median <- function(x){
      median(x[x>0])
    }
    ref_sample <- apply(st_count, 2, geo_mean)
    norm_rle_1 <- sweep(st_count, 2, ref_sample, FUN = "/")
    raw_s_factors <- apply(as.matrix(norm_rle_1), 1, non_zero_median)
    
    ref_sample <- apply(sc_count, 2, geo_mean)
    norm_rle_1 <- sweep(sc_count, 2, ref_sample, FUN = "/")
    raw_v_factors <- apply(as.matrix(norm_rle_1), 1, non_zero_median)
    
    scale_coeff <- exp((-1/(n + m)) * (sum(log(raw_s_factors)) + sum(log(raw_v_factors))))
    scaled_s_factors <- scale_coeff * raw_s_factors
    scaled_v_factors <- scale_coeff * raw_v_factors
  }
  else if(norm_method == "tmm")
  {
    ## TMM(Trimmed Mean Method)
    ### scale_factors
    st_count_t <- t(st_count)
    raw_s_factors <- calcNormFactors(st_count_t, method = "TMM")
    
    sc_count_t <- t(sc_count)
    raw_v_factors <- calcNormFactors(sc_count_t, method = "TMM")
    
    scale_coeff <- exp((-1/(n + m)) * (sum(log(raw_s_factors)) + sum(log(raw_v_factors))))
    scaled_s_factors <- scale_coeff * raw_s_factors
    scaled_v_factors <- scale_coeff * raw_v_factors
  }
  else if(norm_method == "none"){
    scaled_s_factors = rep(1, n)
    scaled_v_factors = rep(1, m)
  }
  else
  {
    stop("Please choose a valid normalization method")
  }
  return(list(s = scaled_s_factors, v = scaled_v_factors))
}


## Get neighbor information ####################################################
vectorized_pdist <- function(A,B){
  an = apply(A, 1, function(rvec) crossprod(rvec,rvec))
  bn = apply(B, 1, function(rvec) crossprod(rvec,rvec))
  m = nrow(A)
  n = nrow(B)
  tmp = matrix(rep(an, n), nrow=m)
  tmp = tmp +  matrix(rep(bn, m), nrow=m, byrow=TRUE)
  sqrt( tmp - 2 * tcrossprod(A,B) )
}

## For 10x Visium data, n_neighbor = 6; For ST data, n_neighbor = 4
get.neighbor <- function(loc, n_neighbor){
  loc <- as.matrix(loc)
  if (n_neighbor == 4){
    loc <- round(loc)
    aa <- sqrt(2)
  } else if (n_neighbor == 6){
    aa <- sqrt(3)
  } else {aa <- 1.2}
  
  dist_matrix <- vectorized_pdist(loc, loc)
  min_dist <- min(dist_matrix[dist_matrix > 0])
  
  dist_threshold <- min_dist*(aa - 1)*0.5 + min_dist
  #print(min_dist)
  #print(dist_threshold)
  
  P <- matrix(0, nrow = nrow(loc), ncol = n_neighbor)
  for (i in 1:nrow(loc)){
    k <- 1
    for (j in 1:nrow(loc)){
      if (dist_matrix[i, j] > 0 & dist_matrix[i, j] < dist_threshold){
        P[i, k] <- j
        k <- k + 1
      }}}
  
  G <- matrix(0, nrow = nrow(loc), ncol = nrow(loc))
  for (i in 1:nrow(loc)) {
    for (j in 1:nrow(loc)) {
      if (dist_matrix[i, j] > 0 & dist_matrix[i, j] < dist_threshold) {
        G[i, j] <- 1
        
      }
    }
  }
  
  return(list("dist" = dist_matrix, "G" = G, "nei" = P))
  
}


## Uses Voronoi tessellation to find nearest neighbours (share a boundary line) and their respective distances.
makeixy = function(data, formula, scale){
  m = model.frame(formula, data=data)
  if(ncol(m)!=3){
    stop("incorrect adjacency formula: id~x+y needed")
  }
  names(m)=c("id","x","y")
  m[,2]=m[,2]/scale
  m[,3]=m[,3]/scale
  m
}

voronoi_adjacency = function(data, formula, scale=1, PLOT=FALSE){
  data = makeixy(data, formula, scale)
  
  P=dim(data)[1];  # number of rows
  
  dd = deldir::deldir(data$x,data$y,suppressMsge=TRUE,plotit=PLOT);  # find adjacencies
  
  ## create adjacency matrix
  A=matrix(0,P,P);
  A[as.matrix(dd$delsgs[,c("ind1","ind2")])] = 1;
  A[as.matrix(dd$delsgs[,c("ind2","ind1")])] = 1;
  
  ## create distance matrix
  D=matrix(NA,P,P);
  D[as.matrix(dd$delsgs[,c("ind1","ind2")])] = sqrt((dd$delsgs[,c("x1")]-dd$delsgs[,c("x2")])^2+(dd$delsgs[,c("y1")]-dd$delsgs[,c("y2")])^2);
  D[as.matrix(dd$delsgs[,c("ind2","ind1")])] = sqrt((dd$delsgs[,c("x1")]-dd$delsgs[,c("x2")])^2+(dd$delsgs[,c("y1")]-dd$delsgs[,c("y2")])^2);
  
  ## create data frame of results
  N=matrix(colSums(A),P,1); # number of adjacencies for each xy$id
  
  ## create neighbor matrix
  n_neighbor <- max(N)
  nei <- matrix(0, nrow = P, ncol = n_neighbor)
  
  for (i in 1:P){
    k <- 1
    for (j in 1:P){
      if (A[i, j] == 1){
        nei[i, k] <- j
        k <- k + 1
      }
    }
  }
  
  return(list(tessellation = dd, dist=D, nei = nei, G = A, NumNeighbours=N, ids=data[,"id"], coords=data[,c("x","y")]));
}


## Visualize metrics results ###################################################
visualize.metrics <- function(data, x, y, fill){
  p = ggplot(data, aes({{x}}, {{y}}, fill={{x}})) +
    geom_boxplot() + xlab("") + ylab("") +
    theme(legend.position="none") +
    facet_grid(~zero) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "grey89", fill=NA, size=0.5))  +
    theme(legend.position="bottom") +
    labs(fill = fill) +
    guides(fill = guide_legend(nrow = 1)) +
    # ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5, size = 10)) 
  
  return(p)
}

visualize.cor <- function(data, x, y, fill){
  p = ggplot(data, aes({{x}}, {{y}}, fill={{x}})) +
    geom_boxplot() + xlab("") + ylab("") +
    theme(legend.position="none") +
    facet_grid(zero~type) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "grey89", fill=NA, size=0.5))  +
    theme(legend.position="bottom") +
    labs(fill = fill) +
    guides(fill = guide_legend(nrow = 1)) +
    # ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5, size = 10))
  print(p)
  
  return(p)
}


## Visualization of clusters ###################################################
plot.domain <- function(data, x, y, size = 4, domain, colors, title = NULL, ncol = 1){
  library(ggplot2)
  p = ggplot(data, aes(x = {{x}}, y = {{y}}, color = {{domain}})) +
    xlab("") + ylab("") +
    geom_point(size=size) +
    coord_fixed(ratio = 1) + 
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())  +
    theme(legend.position="right") +
    theme(plot.margin = margin(0, 0, 0, 0, "cm"),
          panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "grey89", fill=NA, size=0.5)) +
    theme(legend.key=element_blank()) +
    # theme(plot.title = element_text(face="bold")) +
    scale_colour_manual(values = colors) +
    labs(col="Domain") +
    ggtitle(title) +
    guides(color=guide_legend(ncol=ncol))
  
  return(p)
}


## Evaluate results using adjusted Rand index (ARI) ############################
ari <- function(zt, z) {
  library(mclust)
  return(adjustedRandIndex(zt, z))
  
}

## Calculate rmse ##############################################################
get_rmse <- function(actual, predicted){
  rmse = 0
  n <- nrow(actual)
  K <- ncol(actual)
  for (i in 1:n) {
    for (j in 1:K) {
      rmse = rmse + (actual[i, j] - predicted[i, j])^2
    }
  }
  
  rmse = sqrt(rmse/n/K)
  
  return(rmse)
}


## Calculate Matthews correlation coefficient (MCC) ============================
mcc = function(true_gamma, gamma_ppi){
  table = table(true_gamma, gamma_ppi)
  if(length(unique(gamma_ppi)) == 2){
    TN <- table[1]
    FN <- table[2]
    FP <- table[3]
    TP <- table[4]
  }else if (unique(gamma_ppi) == 0){
    TN <- table[1]
    FN <- table[2]
    FP <- 0
    TP <- 0
  }else if (unique(gamma_ppi) == 1){
    TN <- 0
    FN <- 0
    FP <- table[1]
    TP <- table[2]
  } 
  
  TPR = TP/(TP + FN)
  FPR = FP/(FP + TN)
  
  MCC <- (TP * TN - FP * FN)/sqrt((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN))
  
  return(list(table = table, TPR = TPR, FPR = FPR, MCC = MCC))
  
}



## Function to get estimated clusters using pairwise probability matrix (PPM)
get.ppm <- function(z_store, burn, iter, K) {
  n <- ncol(z_store)
  ppm <- matrix(0, nrow = n, ncol = n);
  
  for (ii in (burn + 1):iter) {
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        if (z_store[ii, i] == z_store[ii, j]) {
          ppm[i, j] <- ppm[i, j] + 1;
          ppm[j, i] <- ppm[j, i] + 1;
        }
      }
    }
  }
  
  ppm <- ppm/(iter - burn);
  diag(ppm) <- rep(1, n);
  
  z_ppm <- minbinder(ppm, method = "comp", max.k = K)$cl
  return(list(z_ppm = z_ppm, ppm = ppm))
  
}


## Create the function to get the mode =========================================
getmode <- function(x) {
  uniqx <- unique(x)
  uniqx[which.max(tabulate(match(x, uniqx)))]
}


## Organizing label ============================================================
organize_label = function(z) {
  z_new <- rep(NA, length(z))
  data <- 1
  for(k in unique(z)) {
    z_new[z == k] <- data
    data <- data + 1
  }
  return (z_new)
}


## function to refine the clustering results ###################################
# For an area (e.g., red) with the number of spots <= area_unit, if all neighbors of this area belong to a same cluster (e.g., blue), change the cluster of this area from red to blue. 
cluster_refine <- function(P, cluster, area_unit = 1){
  n_spot <- nrow(P)
  visited <- rep(0, n_spot)
  connected_area <- list()
  k <- 1
  # bfs to find all connected areas
  for (i in 1:n_spot){
    if (visited[i] == 0){
      visited[i] <- 1
      temp <- c(i)
      now_cluster <- cluster[i]
      now_area_list <- c(i)
      while (length(temp) > 0){
        now <- temp[1]
        temp <- temp[-1]
        for (j in P[now, ]){
          if (j != 0){
            if (visited[j] == 0 & cluster[j] == now_cluster){
              visited[j] <- 1
              now_area_list <- c(now_area_list, j)
              temp <- c(temp, j)
            }
          }
        }
      }
      connected_area[[k]] <- now_area_list
      k <- k + 1
    }
  }
  
  n_area <- length(connected_area)
  
  # change the cluster for small areas
  cluster_new <- cluster
  for (i in 1:n_area){
    now_area_list <- connected_area[[i]]
    if (length(now_area_list) <= area_unit){
      # find all neighbors of the current connected area
      neighbor_list <- c()
      for (j in now_area_list){
        neighbor_list <- c(neighbor_list, P[j, P[j, ]!= 0])
      }
      neighbor_list <- setdiff(neighbor_list, now_area_list)
      # cluster of neighbor spots
      neighbor_cluster <- unique(cluster[neighbor_list])
      if (length(neighbor_cluster) == 1){
        cluster_new[now_area_list] <- neighbor_cluster[1]
      }}}
  return(cluster_new)
}


## Calculate Bayesian false discovery rate (BFDR) ##############################
bfdr <- function(PPI, alpha){
  for (c in seq(1,0.1,by=-0.1)) {
    
    BFDR <- sum((1 - PPI)*((1 - PPI) < c))/sum((1 - PPI) < c)
    
    if (BFDR < alpha){
      return(c)
      stop;
    }
  }
}

cal_bfdr <- function(PPI, threshold){
  return(sum((1 - PPI)*((1 - PPI) < threshold))/sum((1 - PPI) < threshold))
}


## Color for scRNA-seq data and ST data ########################################
get_color <- function(data, technology){
  if (data == "STARmap"){
    if (technology == "ST"){
      color <- c("#E64B35FF", "#FFD92F","#4DAF4A", "#7570B3", "#BEAED4", "#FCCDE5", "#377EB8")
    } else {
      color <- c("#CCFF00", "#33FF00", "#FF0000", "#FFD92F", "#00FFFF", "#0066FF", "#CC00FF", "#A6761D", "#BEAED4", "#FCCDE5", "#377EB8")
    }
  } else if (data == "MOB"){
    if (technology == "ST"){
      color <- c("#E64B35FF", "#FFD92F","#4DAF4A", "#7570B3")
    } else {
      color <- c("#BEAED4", "#FCCDE5", "#377EB8","#FDC086", "#7FC97F")
    }
  }
  
  return(color)
}


## Function to plot trace figure and PPI
plot.trace <- function(data, x, y, ylab) {
  p <- ggplot(data, aes(x={{x}}, y={{y}})) + 
    geom_line() +
    ylab(ylab) +
    theme(legend.key=element_rect(fill="white")) +
    theme_minimal() +
    theme(panel.grid = element_blank(), 
          panel.border = element_rect(fill= "transparent")) +
    theme(legend.key.size = unit(0.5, "cm"))
  
  return(p)
}

plot.ppi <- function(data_type, data, x, y, z=NULL) {
  if (data_type == "real") {
    p <- ggplot(data, aes(x={{x}},xend={{x}},y=0, yend={{y}})) + 
      geom_segment()
      
  } else if (data_type == "simulated"){
    p <- ggplot(data, aes(x={{x}}, y={{z}})) + 
      geom_point(color="red") +
      geom_segment(aes(x={{x}},xend={{x}},y=0, yend={{y}}))
  }
  
  p <- p + theme(legend.key=element_rect(fill="white")) +
    theme_minimal() +
    theme(panel.grid = element_blank(), 
          panel.border = element_rect(fill= "transparent")) +
    geom_hline(yintercept = 0.5, colour = "red", linetype=3) +
    geom_hline(yintercept = 1-c, colour = "blue", linetype=3) +
    xlab("Gene index") +
    ylab("Posterior probabilities of inclusion") + 
    theme(plot.caption = element_text(hjust=0.5))
  
  return(p)
  
}


