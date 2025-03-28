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

library(reshape2)
library(gtools)
library(scatterpie)
library(Seurat)
library(ggpubr)

## Data quality control (QC) ###################################################
## scRNA-seq QC
sc.qc <- function(count, metadata) {
  ## Sample-wise quality control
  index <- which(rowSums(count > 0) > 0)
  count <- count[index, ]
  metadata <- metadata[index, ]
  
  ## Feature-wise quality control
  index <- which(colSums(count > 0) > 0)
  count <- count[, index]
  
  return (list("count" = count, "metadata" = metadata))
}

## ST QC
st.qc <- function(count, loc, cutoff_sample = 100, cutoff_feature = 0.1) {
  ## Sample-wise quality control
  index <- which(rowSums(count) >= cutoff_sample)
  count <- count[index, ]
  loc <- loc[index, ]
  
  ## Feature-wise quality control
  index <- which(colSums(count > 0) >= cutoff_feature*nrow(loc))
  if (sum(count==0)/nrow(count)/ncol(count) > 0.9) {
    index <- which(colSums(count > 0) > 0)
  }
  count <- count[, index]
  
  return (list("count" = count, "loc" = loc))
}


## Data normalization #######################################################
normalize.separate <- function(count, norm_method = 'tss'){
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


## Find marker genes ###########################################################
find.markers <- function(count, metadata, min.pct=0.1, logfc.threshold=0.25, test.use = "wilcox"){
  # Create Seurat object from raw counts matrix
  seu <- CreateSeuratObject(counts = t(as(as.matrix(count), "dgCMatrix")))
  
  # Add clustering information to Seurat object
  seu <- AddMetaData(seu, metadata = metadata, col.name = "cluster")
  
  # Set the identity classes (clusters)
  Idents(seu) <- seu$cluster
  
  # Normalize the data
  seu <- NormalizeData(seu)
  
  # Identify the 2000 most variable features
  seu <- FindVariableFeatures(seu, selection.method = "vst")
  
  # Scale the data
  seu <- ScaleData(seu)
  
  # Run PCA
  seu <- RunPCA(seu, features = VariableFeatures(object = seu))
  
  # Find markers differential each group
  seu.markers <- FindAllMarkers(seu, min.pct = min.pct, logfc.threshold = logfc.threshold, only.pos=TRUE,test.use = test.use)
  seu.markers <- seu.markers[seu.markers$p_val_adj<0.05, ]
  
  return(seu.markers)
}


## Visualization of domains ###################################################
plot.domain <- function(loc, size = 4, domain, label = "Domain", colors=NULL, title = NULL, ncol = 1){
  library(ggplot2)
  data=data.frame("x"=loc[, 1],"y"=loc[, 2],"domain"=as.factor(domain))
  
  if(is.null(colors)){
    colors = rainbow(length(unique(data$domain)))
  }else{
    colors = colors
  }
  p = ggplot(data, aes(x = x, y = y, color = domain)) +
    labs(col=label, title = title, xlab = "", ylab = "") +
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
          panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
    theme(legend.key=element_blank()) +
    # theme(plot.title = element_text(face="bold")) +
    scale_colour_manual(values = colors) +
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

get_pcc <- function(actual, predicted){
  pcc = cor(c(as.matrix(actual)), c(as.matrix(predicted)))
  return(pcc)
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



## Function to plot trace figure and PPI
plot.trace <- function(x, y, ylab, ylim) {
  data <- data.frame("x" = x, "y" = y)
  p <- ggplot(data, aes(x={{x}}, y={{y}})) + 
    geom_line() +
    xlab("Iteration") +
    ylab(ylab) +
    theme(legend.key=element_rect(fill="white")) +
    theme_minimal() +
    theme(panel.grid = element_blank(), 
          panel.border = element_rect(fill= "transparent")) +
    theme(legend.key.size = unit(0.5, "cm")) +
    ylim(ylim)
  
  return(p)
}

plot.trace.chain <- function(data) {
  ggplot(data, aes(x = Iteration, y = Value, group = Chain, color = Chain)) +
    theme_bw() +
    geom_line() +
    theme(panel.grid = element_blank()) +
    scale_color_manual(values=c("#1B9E77", "#D95F02", "#7570B3", "#E7298A")) +
    theme(legend.position = c(0.995, 0.003),
          legend.justification = c("right", "bottom"),
          legend.key.size = unit(1, "lines")) +
    ylab("1st Spot largest cell type proportion") 
}


plot_heatmap <- function(x, y, fill, xlab=NULL, ylab=NULL){
  data <- data.frame("x"= x , "y" = y)
  ggplot(data, aes(x=x, y=y, fill = {{fill}})) +
    geom_tile()+
    scale_fill_gradient2(low = "blue", 
                         mid = "white",
                         high = "red",  
                         midpoint = 0, 
                         limit = c(-1,1), 
                         space = "Lab", 
                         name="Pearson\ncorrelation") +
    theme(plot.margin = unit(c(0,0,0,0), "cm")) +
    # theme_minimal()+ 
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white"))+
    coord_fixed() +
    xlab(xlab) +
    ylab(ylab) 
}


plot.prop <- function(proportion,spatial_location,ct.visualize = ct.visualize,colors = c("lightblue","lightyellow","red"),NumCols, pointSize = 3.0){
  if(is.null(colors)){
    colors = c("lightblue","lightyellow","red")
  }else{
    colors = colors
  }
  res_CARD = as.data.frame(proportion)
  res_CARD = res_CARD[,order(colnames(res_CARD))]
  location = as.data.frame(spatial_location)
  if(sum(rownames(res_CARD)==rownames(location))!= nrow(res_CARD)){
    stop("The rownames of proportion data does not match with the rownames of spatial location data")
  }
  ct.select = ct.visualize
  res_CARD = res_CARD[,ct.select]
  if(!is.null(ncol(res_CARD))){
    res_CARD_scale = as.data.frame(apply(res_CARD,2,function(x){
      (x - min(x)) / (max(x) - min(x))
    } ))}else{
      res_CARD_scale = as.data.frame((res_CARD - min(res_CARD)) / (max(res_CARD) - min(res_CARD)))
      colnames(res_CARD_scale) = ct.visualize
    }
  res_CARD_scale$x = as.numeric(location$x)
  res_CARD_scale$y = as.numeric(location$y)

  mData = melt(res_CARD_scale,id.vars = c("x","y"))
  colnames(mData)[3] <- "Cell_Type"
  b = c(0,1)
  
  p = suppressMessages(ggplot(mData, aes(x, y)) + 
                         geom_point(aes(colour = value),size = pointSize) +
                         scale_color_gradientn(colours = colors) + 
                         #scale_color_viridis_c(option = 2)+
                         scale_x_discrete(expand = c(0, 1)) + scale_y_discrete(expand = c(0,1))+ 
                         facet_wrap(~Cell_Type, ncol = NumCols)+
                         # facet_grid(Method~Cell_Type)+  
                         coord_fixed()+
                         theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
                               #legend.position=c(0.14,0.76),
                               panel.background = element_blank(),
                               plot.background = element_blank(),
                               panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
                               axis.text =element_blank(),
                               axis.ticks =element_blank(),
                               axis.title =element_blank(),
                               legend.title=element_text(size = 10),
                               legend.text=element_text(size = 10),
                               strip.text = element_text(size = 10),
                               legend.key = element_rect(colour = "transparent", fill = "white"),
                               legend.key.size = unit(0.45, 'cm'))) +
    labs(color = "Proportion")
  return(list("p" = p, "data" = mData))
}


plot.pie <- function(proportion,spatial_location,colors = NULL,radius = NULL,seed = NULL, nrow=2){
  res_CARD = as.data.frame(proportion)
  res_CARD = res_CARD[,mixedsort(colnames(res_CARD))]
  location = as.data.frame(spatial_location)
  if(sum(rownames(res_CARD)==rownames(location))!= nrow(res_CARD)){
    stop("The rownames of proportion data does not match with the rownames of spatial location data")
  }
  colorCandidate = c("#1e77b4","#ff7d0b","#ceaaa3","#2c9f2c","#babc22","#d52828","#9267bc",
                     "#8b544c","#e277c1","#d42728","#adc6e8","#97df89","#fe9795","#4381bd","#f2941f","#5aa43a","#cc4d2e","#9f83c8","#91675a",
                     "#da8ec8","#929292","#c3c237","#b4e0ea","#bacceb","#f7c685",
                     "#dcf0d0","#f4a99f","#c8bad8",
                     "#F56867", "#FEB915", "#C798EE", "#59BE86", "#7495D3",
                     "#D1D1D1", "#6D1A9C", "#15821E", "#3A84E6", "#997273",
                     "#787878", "#DB4C6C", "#9E7A7A", "#554236", "#AF5F3C",
                     "#93796C", "#F9BD3F", "#DAB370", "#877F6C", "#268785",
                     "#f4f1de","#e07a5f","#3d405b","#81b29a","#f2cc8f","#a8dadc","#f1faee","#f08080")
  if(is.null(colors)){
    #colors = brewer.pal(11, "Spectral")
    if(ncol(res_CARD) > length(colorCandidate)){
      colors = colorRampPalette(colorCandidate)(ncol(res_CARD))
    }else{
      
      if(is.null(seed)){
        iseed = 12345
      }else{
        iseed = seed
      }
      set.seed(iseed)
      colors = colorCandidate[sample(1:length(colorCandidate),ncol(res_CARD))]
    }
  }else{
    colors = colors
  }
  data = cbind(res_CARD,location)
  ct.select = colnames(res_CARD)
  
  if(is.null(radius)){
    radius = (max(data$x) - min(data$x)) * (max(data$y) - min(data$y))
    radius = radius / nrow(data)
    radius = radius / pi
    radius = sqrt(radius) * 0.85
  }else{
    #### avoid the situation when the radius does not generate the correct figure
    radius = radius
  }
  p = suppressMessages(ggplot() + geom_scatterpie(aes(x=x, y=y,r = radius),data=data,
                                                  cols=ct.select,color=NA) + 
                         # coord_fixed(ratio = 1*max(data$x)/max(data$y)) + 
                         coord_fixed(ratio = 1) + 
                         scale_fill_manual(values =  colors)+
                         theme(plot.margin = margin(0, 0, 0, 0, "cm"),
                               panel.background = element_blank(),
                               plot.background = element_blank(),
                               panel.border = element_rect(colour = "black", fill=NA, size=0.5),
                               axis.text =element_blank(),
                               axis.ticks =element_blank(),
                               axis.title =element_blank(),
                               legend.title=element_text(size = 16),
                               legend.text=element_text(size = 15),
                               legend.key = element_rect(colour = "transparent", fill = "white"),
                               # legend.key.size = unit(0.45, 'cm'),
                               strip.text = element_text(size = 16),
                               legend.position="right")+
                         guides(fill=guide_legend(title="Cell type", nrow = nrow))
                         ) 
  return(p)
}


## Visualization correlation coefficient matrix ################################
plot.cor <- function(proportion) {
  cor_mat <- cor(as.matrix(proportion))
  
  data=expand.grid(X=colnames(cor_mat),Y=colnames(cor_mat))
  data$Z <- expand.grid(cor_mat)$Var1
  
  p = ggplot(data, aes(factor(X), factor(Y), fill= Z)) +  
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", 
                         high = "red", 
                         mid = "white", 
                         midpoint = 0, 
                         limit = c(-1, 1),
                         space = "Lab",
                         na.value = "red",
                         name="Correlation") +
    theme_bw()+
    labs(y=" ",x=" ",title =" ")+
    theme(axis.ticks.x  = element_blank(), 
          axis.ticks.y  = element_blank(), 
          axis.text.x= element_text(size = 16,hjust = 1,angle = 50,color="black"), 
          axis.text.y = element_text(size = 16,color="black"),
          legend.text = element_text(size = 16,angle= 0.5,color="black"), 
          legend.title = element_text(size = 16,color="black"), #element_blank(),  #element_text(size = 14),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_blank())    +
    theme(plot.margin = unit(c(0,0,0,0), "cm"))

  
  return(p)
  
}


## Get initial of clustering from Louvain method ###############################
cluster_start <- function(count, N_domain) {
  if (is.null(colnames(count))) {
    colnames(count) = paste0("gene", 1:p)
  }
  if (is.null(rownames(count))) {
    rownames(count) = paste0("sample", 1:n)
  }
  
  colData <- data.frame("sample" = rownames(count))
  rowData <- data.frame("gene" = colnames(count))
  
  sce <- SingleCellExperiment(assays = list(counts = as(t(count), 'dgCMatrix')), rowData = rowData, colData = colData)
  log_count <- log(sce@assays@data@listData$counts + 1) 
  suppressMessages(M1 <- as(log_count, "dgCMatrix"))
  
  sce@assays@data@listData$logcounts <- M1
  rm(log_count)
  
  ss <- as.Seurat(sce, counts = "counts", data = "logcounts", assay = NULL, project = "SingleCellExperiment")
  ss <- SCTransform(ss, assay = "originalexp", verbose = FALSE)
  ss <- RunPCA(ss, assay = "SCT", verbose = FALSE)
  ss <- FindNeighbors(ss, verbose = FALSE)
  
  i <- 0.1
  
  while (i <= 1 & i >= 0.1) {
    ss2 <- FindClusters(ss, verbose = FALSE, resolution = i)
    result <- ss2@meta.data
    if (length(unique(result$seurat_clusters)) == N_domain){
      break
    }
    i = i + 0.1
  }
  
  z <- as.numeric(result$seurat_clusters)-1
  
  if (length(unique(result$seurat_clusters)) != N_domain){
    set.seed(12345)
    z = kmeans(count, centers = K)$cluster
    z = z - 1
  }
  return(z)
}


## Suppressing output from cat() ###############################################
quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 
