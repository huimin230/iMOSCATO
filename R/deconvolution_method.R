################################################################################
## 
## Title        : Robust Bayesian Integrative Modeling of Single Cell and Spatially Resolved Transcriptomics Data
## Project      : iMOSCATO
## Research goal: Bayesian cell-type deconvolution with spatial domain detection
## Authors      : Huimin Li
## Contact      : huimin.li@utdallas.edu
## Code         : Cell-type deconvolution methods
## Modified     : 2025-04-01
## 
################################################################################


## Run CARD ####################################################################
## Tutorial: https://yingma0107.github.io/CARD/documentation/04_CARD_Example.html
CARD_run <- function(sc_count, sc_meta, st_count, loc, minCountGene = 0, minCountSpot = 0) {
  library(CARD)
  cat("Running CARD ....\n")
  
  start_time = Sys.time()
  cellType <- sort(unique(sc_meta$cellType))
  
  if(!("sampleInfo"%in%names(sc_meta))) {
    sc_meta$sampleInfo <- rep("sample1", nrow(sc_meta))
  }
  
  ## Create CARD object
  CARD_object <- createCARDObject(
    sc_count = t(as.matrix(sc_count)),
    sc_meta = sc_meta,
    spatial_count = t(as.matrix(st_count)),
    spatial_location = loc,
    ct.varname = "cellType",
    ct.select = cellType,
    sample.varname = "sampleInfo",
    minCountGene = 100,
    minCountSpot = 5)
  
  ## Run CARD
  CARD_object <- suppressWarnings(CARD_deconvolution(CARD_object))
  
  ## Save results
  prop_result <- as.data.frame(CARD_object@Proportion_CARD)
  rownames(prop_result) <- rownames(CARD_object@spatial_location)
  prop_result <- prop_result[, cellType]
  
  B = t(CARD_object@algorithm_matrix$B)
  marker_gene_result <- data.frame("gene"= rownames(CARD_object@algorithm_matrix$B))
  
  end_time = Sys.time()
  time <- difftime(end_time, start_time, units = "secs")
  
  cat("Running time for CARD:", round(time, digits = 2), "seconds","\n")
  
  return(list("prop_result"=prop_result, "reference_matrix" = B,
              "marker_gene_result"=marker_gene_result, "time"=time))
}


## Run RCTD ####################################################################
## Tutorial: https://raw.githack.com/dmcable/spacexr/master/vignettes/spatial-transcriptomics.html
RCTD_run <- function(sc_count, sc_meta, st_count, loc, CELL_MIN_INSTANCE=25,min_UMI = 100, UMI_min=100, require_int=TRUE) {
  library(spacexr)
  cat("Running RCTD ....\n")
  
  start_time = Sys.time()

  ## Single-Cell Reference
  counts <- t(as.matrix(sc_count))
  cellType <- sort(unique(sc_meta$cellType))
  
  ## Change cell type name if contains "/"
  cell_types <- sc_meta$cellType
  cell_types <- gsub("/", "-", cell_types)
  cell_types <- as.factor(as.character(cell_types))
  names(cell_types) <- rownames(sc_meta)
  
  if ("nUMI"%in%names(sc_meta)) {
    nUMI <- sc_meta$nUMI
  } else {
    nUMI <- colSums(counts)   # In this case, total counts is nUMI
  }
  
  names(nUMI) <- rownames(sc_meta)
  
  ## Create the Reference object
  reference <- Reference(counts=counts, 
                         cell_types=cell_types, 
                         nUMI=nUMI,
                         min_UMI = min_UMI,
                         require_int = require_int)

  ## Create SpatialRNA object
  counts <- t(as.matrix(st_count))
  coords <- loc
  nUMI <- colSums(counts) # In this case, total counts per pixel is nUMI
  puck <- SpatialRNA(coords=coords, 
                              counts=counts, 
                              nUMI=nUMI)
  
  # ## Create RCTD object
  # myRCTD <- suppressMessages(create.RCTD(spatialRNA=puck, 
  #                                                 reference=reference, 
  #                                                CELL_MIN_INSTANCE = CELL_MIN_INSTANCE))
  
  myRCTD <- create.RCTD(spatialRNA=puck, 
                        reference=reference, 
                        CELL_MIN_INSTANCE = CELL_MIN_INSTANCE,
                        UMI_min = UMI_min)
  
  
  ## Run RCTD
  myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')
  
  ## Save results
  results <- myRCTD@results
  
  norm_weights <- normalize_weights(results$weights) 
  prop_result <- as.data.frame(as.matrix(norm_weights))
  
  if (length(grep("/", cellType))) {
    names(prop_result) <- gsub("-", "/", names(prop_result))
  }
  prop_result <- prop_result[, cellType]
  
  marker_gene_result <- data.frame("gene" = myRCTD@internal_vars$gene_list_bulk)
  
  end_time = Sys.time()
  time <- difftime(end_time, start_time, units = "secs")
  
  cat("Running time for RCTD:", round(time, digits = 2), "seconds","\n")
  
  return(list("prop_result"=prop_result, "marker_gene_result"=marker_gene_result, "time"=time))
}


# ## Run SpatialDWS ##############################################################
# ## Tutorial: https://rubd.github.io/Giotto_site/articles/tut7_giotto_enrichment.html#hypergeometric-enrichment
# SpatialDWS_run <- function() {
#   library(Giotto)
#   cat("Running SpatialDWS ....\n")
#   
#   start_time = Sys.time()
#   cellType <- sort(unique(sc_meta$cellType))
#   
#   ## 1. processing steps
#   path_to_matrix = t(st_count) 
#   path_to_locations = loc
#   my_giotto_object = createGiottoObject(raw_exprs = path_to_matrix,
#                                         spatial_locs = path_to_locations)
#   
#   ## processing
#   my_giotto_object <- normalizeGiotto(gobject = my_giotto_object)
#   
#   ## 2. Run spatial cell type enrichments methods
#   ## 2.1 PAGE enrichment
#   ## highly variable genes (HVG)
#   count = sc_count
#   metadata = sc_meta[, "cellType"]
#   
#   full.seu.markers <- find.markers(count=count,
#                                    metadata=metadata,
#                                    min.pct = 0.1,
#                                    logfc.threshold = 1.25,
#                                    test.use = "wilcox")
#   
#   seu.markers = full.seu.markers
#   table(seu.markers$cluster)[cellType]
#   gene_all <- list()
#   
#   for (k in 1:length(cellType)) {
#     data <- seu.markers[seu.markers$p_val_adj<0.05 & seu.markers$cluster==cellType[k], ]
#     data <- data[order(data$p_val_adj), ]
#     gene_all[[k]] <- data$gene
#   }
#   
#   signature_matrix = makeSignMatrixPAGE(sign_names = cellType, sign_list = gene_all)
#   
#   # runSpatialEnrich() can also be used as a wrapper for all currently provided enrichment options
#   my_giotto_object = runPAGEEnrich(gobject = my_giotto_object,
#                                    sign_matrix = signature_matrix,
#                                    min_overlap_genes = 2)
#   
#   ## 2.2 RANK enrichment
#   single_cell_matrix = t(sc_count)
#   cell_annotations = sc_meta$cellType
#   
#   ## 1.2 create rank matrix
#   rank_matrix = makeSignMatrixRank(sc_matrix = single_cell_matrix, sc_cluster_ids = cell_annotations)
#   
#   ## 1.3 enrichment test with RANK runSpatialEnrich() can also be used as a wrapper for all currently provided enrichment options
#   my_giotto_object = runRankEnrich(gobject = my_giotto_object, sign_matrix = rank_matrix)
#   
#   ## 2.3 Hypergeometric enrichment
#   my_giotto_object = runHyperGeometricEnrich(gobject = my_giotto_object,
#                                              sign_matrix = signature_matrix)
#   
#   cell_types_subset = colnames(signature_matrix)
#   spatCellPlot(gobject = my_giotto_object, 
#                spat_enr_names = 'hypergeometric',
#                cell_annotation_values = cell_types_subset,
#                cow_n_col = 2,coord_fix_ratio = NULL, point_size = 2.75)
#   
#   ## 2.4 Deconvolution
#   my_giotto_object@cellType = sc_meta$cellType
#   my_giotto_object = runDWLSDeconv(gobject = my_giotto_object, 
#                                    cluster_column = "cellType",
#                                    sign_matrix = signature_matrix)
#   
#   prop_result <- prop_result[, cellType]
#   
#   end_time = Sys.time()
#   time <- difftime(end_time, start_time, units = "secs")
#   
#   cat("Running time for SpatialDWS:", round(time, digits = 2), "seconds","\n")
#   
#   return(list("prop_result"=prop_result, "marker_gene_result"=marker_gene_result, "time"=time))
#   
# }
# 
# 
# 
# 
# ## Run SPOTlight ##############################################################
# SPOTlight_run <- function() {
#   library(CARD)
#   cat("Running SpatialDWS ....\n")
#   
#   start_time = Sys.time()
#   cellType <- sort(unique(sc_meta$cellType))
#   
#   prop_result <- prop_result[, cellType]
#   
#   end_time = Sys.time()
#   time <- difftime(end_time, start_time, units = "secs")
#   
#   
#   cat("Running time for SPOTlight:", round(time, digits = 2), "seconds","\n")
#   
#   return(list("prop_result"=prop_result, "marker_gene_result"=marker_gene_result, "time"=time))  
# }