################################################################################
## 
## Paper        : Bayesian Integrative Modeling of Single Cell and Spatial Transcriptomics Data
## Project      : iMOSCATO 
## Research goal: Bayesian cell type deconvolution with feature selection
## Modified     : 2023-07-12
## 
################################################################################

################################################################################
## 
## Paper 1      : Spatially Informed Cell Type Deconvolution for Spatial Transcriptomics
## Method       : CARD
## Tutorial     : https://yingma0107.github.io/CARD/documentation/04_CARD_Example.html
## 
## Paper 2      : Robust decomposition of cell type mixtures in spatial transcriptomics
## Method       : RCTD
## Tutorial     : https://raw.githack.com/dmcable/spacexr/master/vignettes/spatial-transcriptomics.html
##
################################################################################


rm(list = ls())


## Write functions to implement cell type deconvolution methods ################
deconvolution <- function(method, X, Y, loc, sc_meta, st_meta, true_prop = NULL){
  library(CARD)  
  library(spacexr)
  library(Matrix)
  library(ggplot2)
  
  sc_meta$cellType <- factor(sc_meta$cellType)
  cellType <- levels(sc_meta$cellType)
  
  if (length(true_prop) > 0){
    a <- which(names(true_prop) == cellType[1])
  }
  
  
  if (method == "CARD") {
    ## Create CARD object ######################################################
    sc_count = t(as.matrix(X))
    st_count = t(as.matrix(Y))
    ct.select = cellType
    start_time = proc.time()
    
    CARD_object = createCARDObject(
      sc_count = sc_count,
      sc_meta = sc_meta,
      spatial_count = st_count,
      spatial_location = loc,
      ct.varname = "cellType",
      ct.select = ct.select,
      sample.varname = "sampleInfo")
    
    
    ## Running CARD ############################################################
    CARD_object = CARD_deconvolution(CARD_object)
    
    end_time = proc.time()
    time = as.numeric((end_time - start_time)[3], "secs")
    
    
    ## CARD results ############################################################
    print(round(apply(CARD_object@Proportion_CARD, 2, summary),digits = 4))
    
    ## Estimated cell type proportions 
    p_hat <- as.data.frame(CARD_object@Proportion_CARD)
    prop_data <- p_hat
    names(prop_data) <- cellType
    prop_data$x <- CARD_object@spatial_location$x
    prop_data$y <- CARD_object@spatial_location$y
    
    results <- CARD_object
    marker_gene_data <- data.frame("gene"= rownames(CARD_object@algorithm_matrix$B))

  } else if (method == "RCTD") {
    ## Single-Cell Reference ###################################################
    sc_count <- X
    counts <- t(as.matrix(sc_count))
    cell_types <- as.character(sc_meta$cellType)
    names(cell_types) <- rownames(sc_meta)
    cell_types <- gsub("/", "-", cell_types)
    table(cell_types)
    cell_types <- as.factor(cell_types)
    
    if ("nUMI"%in%names(sc_meta)) {
      nUMI <- sc_meta$nUMI
    } else {
      nUMI <- colSums(counts)   # In this case, total counts is nUMI
    }
    
    names(nUMI) <- rownames(sc_meta)
    
    ## Create the Reference object
    reference <- Reference(counts, cell_types, nUMI, require_int = FALSE)
    
    ## Examine reference object (optional)
    print(dim(reference@counts)) #observe Digital Gene Expression matrix
    table(reference@cell_types) #number of occurences for each cell type
    
    
    ## Spatial Transcriptomics data ############################################
    st_count <- Y
    counts <- t(as.matrix(st_count))
    coords <- loc
    nUMI <- colSums(counts) # In this case, total counts per pixel is nUMI
    
    ## Create SpatialRNA object
    puck <- SpatialRNA(coords, counts, nUMI)
    
    ## Examine SpatialRNA object (optional)
    print(dim(puck@counts)) # observe Digital Gene Expression matrix
    hist(log(puck@nUMI,2)) # histogram of log_2 nUMI
    
    print(head(puck@coords)) # start of coordinate data.frame
    barcodes <- colnames(puck@counts) # pixels to be used (a list of barcode names). 
    
    # This list can be restricted if you want to crop the puck e.g. 
    # puck <- restrict_puck(puck, barcodes) provides a basic plot of the nUMI of each pixel
    # on the plot:
    plot_puck_continuous(puck, barcodes, puck@nUMI, 
                         ylimit = c(0,round(quantile(puck@nUMI,0.9))), 
                         title ='plot of nUMI') 

    
    ## Creating RCTD Object ####################################################
    myRCTD <- create.RCTD(puck, reference)
    marker_gene_data <- data.frame("gene" = myRCTD@internal_vars$gene_list_bulk)
    
    
    ## Running RCTD ############################################################
    start_time = proc.time()
    myRCTD <- run.RCTD(myRCTD, doublet_mode = "doublet")
    
    end_time = proc.time()
    time = as.numeric((end_time - start_time)[3], "secs")
    
    
    ## RCTD results ############################################################
    results <- myRCTD@results
    
    ## Normalize the cell type proportions to sum to 1.
    norm_weights = normalize_weights(results$weights) 
    rowSums(norm_weights)
    p_hat = norm_weights
    p_hat <- as.matrix(p_hat)
    
    # cell_type_names <- myRCTD@cell_type_info$info[[2]] #list of cell type names
    # spatialRNA <- myRCTD@spatialRNA
    # resultsdir <- 'RCTD_Plots' ## you may change this to a more accessible directory on your computer.
    # dir.create(resultsdir)
    
    ## Estimated cell type proportions #########################################
    prop_data <- as.data.frame(p_hat)
    names(prop_data) <- cellType
    prop_data$x <- loc$x
    prop_data$y <- loc$y
    
  }
  
  ## Estimated dominant cell type ############################################
  prop_data$dominant_hat <- sapply(1:dim(p_hat)[1], function(x) colnames(prop_data[, cellType])[which.max(prop_data[, cellType][x, ])])
  
  prop_data <- merge(st_meta, prop_data, by = c("x", "y"))
  ARI <- ari(prop_data$Layer, prop_data$dominant_hat)
  
  scatterpie <- CARD.visualize.pie(proportion = prop_data[, cellType] , spatial_location = prop_data[, c("x", "y")], colors = sc_color) +
    guides(fill=guide_legend(nrow=2, byrow = TRUE)) +                         
    labs(fill="Cell type") +
    theme(strip.text = element_text(size = 10,face="bold"))
  
  
  ## Save result ###############################################################
  if (length(true_prop) > 0) {
    temp <- merge(true_prop, prop_data, by = names(true_prop)[1:(a-1)])
    actual <- temp[, paste0(cellType, ".x")]
    predicted <- temp[, paste0(cellType, ".y")]
    RMSE <- get_rmse(actual = actual, predicted = predicted)
    PCC <- cor(c(as.matrix(actual)), c(as.matrix(predicted)))
    
    output = data.frame("data" = data_name, "method" = method, 
                        "metrics" = c("RMSE", "PCC", "ARI", "time"), 
                        "value" = c(RMSE, PCC, ARI, time))
  } else {
    output = data.frame("data" = data_name, "method" = method, 
                        "metrics" = c("ARI", "time"), 
                        "value" = c(ARI, time))
    
  }

  
  return(list("results"=results, "prop_data"=prop_data, "marker_gene_data"=marker_gene_data, "output"=output, "scatterpie"=scatterpie))
  
}
