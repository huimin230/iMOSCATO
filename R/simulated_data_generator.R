################################################################################
## 
## Title        : Robust Bayesian Integrative Modeling of Single Cell and Spatially Resolved Transcriptomics Data
## Project      : iMOSCATO
## Research goal: Bayesian cell-type deconvolution with spatial domain detection
## Authors      : Huimin Li
## Contact      : huimin.li@utdallas.edu
## Code         : Generate simulated data
## Modified     : 2025-04-01
## 
################################################################################

rm(list = ls())

## Write function to generate scRNA-seq data ===================================
simulate.data = function(M=500, P = 1000, P_gamma = 20, K, N_k, zero=0.05, sparsity=0.25, seed = NA){
  if (!is.na(seed)) {
    set.seed(seed)
  }
  
  cat(paste0("Simulating cell-level data! \n"))
  
  ## size factor
  v = exp(rnorm(M, mean = 0, sd = 0.2))
  
  ## dispersion parameter
  psi = rexp(P, rate = 1/psi_mean)

  ## genes
  gamma = rep(FALSE, P)
  gamma[sample(1:P, P_gamma)] = TRUE
  
  ## cell name
  cell_name = unlist(sapply(1:K, function(k) rep(paste0("cellType", k), N_k[k])))
  table(cell_name)
  
  ## random errors
  sigma = rexp(P, rate = 1/sigma_mean)
  varepsilon = matrix(rnorm(P, mean = 0, sd = sigma), nrow = M, ncol = P, byrow = TRUE)
  
  ## log normalized expression levels with effect size
  loglambda = beta_0 + varepsilon
  temp = loglambda
  
  for (j in 1:P_gamma) {
    e_k_temp = (e_k[sample(1:K, K)])
    beta_1 = unlist(sapply(1:K, function(k) rep(e_k_temp[k], N_k[k])))
    table(beta_1)
    loglambda[, gamma][, j] = loglambda[, gamma][, j] + beta_1
  }
  
  
  ## generate count matrix
  count = matrix(rnbinom(M * P,
                         mu = v * exp(loglambda), 
                         size = matrix(psi, nrow = M, ncol = P, byrow = TRUE)), 
                 nrow = M, ncol = P)
  
  ## zero inflation settings
  Zeta = matrix(0, nrow = M, ncol = P)
  
  if (zero > 0) {
    Zeta = matrix(rbinom(M * P, 1, zero), nrow = M, ncol = P)
    count[which(Zeta == 1)] = 0
  }
  
  ## normalized level matrix: K-by-P not N-by-P to reduce dimension
  sc_B = exp(loglambda)
  B = sc_B[!duplicated(sc_B), ]
  rownames(B) = paste0("cellType", 1:K)

  sc_meta = data.frame("Sample" = 1:nrow(count), "cellType" = cell_name)
  rownames(sc_meta) <- paste0("Sample", 1:nrow(count))
  
  rownames(sc_B)=rownames(count)=rownames(Zeta)=rownames(sc_meta)
  
  sc_count <- count
  
  ## generate group-specific cell type proportion matrix
  cat(paste0("Simulating cell type proportions! \n"))
  
  ## number of locations
  N = nrow(loc)
  
  ## spatial domain allocation parameter
  rownames(loc) <- paste0(round(loc$x), "x", round(loc$y))
  sum(duplicated(rownames(loc)))
  
  ## zero-inflation indicator of Omega
  Delta = matrix(0, nrow = N, ncol = K)
  
  for (domain in 1:2) {
    if (domain == 1){
      dominant = 1
    }else{
      dominant = 4
    }
    
    Delta[which(st_meta$Domain==domain), -dominant] = matrix(rbinom(sum(st_meta$Domain==domain)*(K-1), 1, sparsity), nrow = sum(st_meta$Domain==domain), ncol = K-1)
  }
  
  
  ## concenration parameters
  A <- matrix(c(3, 1, 1, 1, 1, 1, 1, 3), byrow = TRUE, nrow = 2)
  
  ## create empty cell type proportion matrix
  Omega = matrix(NA, nrow = N, ncol = K)
  
  for (i in 1:N) {
    for (k in 1:K) {
      if (Delta[i, k] == 1) {## if Delta_{ik} = 0, Pi_{ik} = Omega_{ik} = 0
        Omega[i, k] = 0
        
      } else{
        Omega[i, k] = rgamma(1, A[st_meta$Domain[i], k], 1)
        
      }
      
    } ## end of k
  } ## end of i
  
  ## To prevent all cell types have zero proportion in a spot
  index <- which(rowSums(Omega)==0)
  
  if (length(index) > 0) {
    k <- sapply(index, function(x) which.max(A[st_meta$Domain[x], ]))
    Omega[index, k] <- 1
    
  }
  
  if (length(which(rowSums(Omega)==0)) > 0) {
    stop("Cell type proportions cannot all be zero in a spot!")
  }
  
  Pi = Omega/rowSums(Omega)
  
  colnames(Omega) = colnames(Pi) = colnames(Delta) = colnames(A) = ct.select
  
  cat(paste0("Simulating spot-level data! \n"))
  
  ## Set parameters
  ## size factor
  s = exp(rnorm(N, mean = 0, sd = 0.2))
  
  ## dispersion parameter
  phi = rexp(P, rate = 1/phi_mean)
  
  ## calcuate weighted average normalized gene expressions
  st_B = as.matrix(Pi)%*%as.matrix(B)
  
  ## generate count matrix
  count = matrix(rnbinom(N * P,
                         mu = s * st_B, 
                         size = matrix(phi, nrow = N, ncol = P, byrow = TRUE)), 
                 nrow = N, ncol = P)
  
  ## zero inflation settings
  Eta = matrix(0, nrow = N, ncol = P)
  
  if (zero > 0) {
    Eta = matrix(rbinom(N * P, 1, zero), nrow = N, ncol = P)
    count[which(Eta == 1)] = 0
  }
  
  st_count <- count
  
  rownames(st_meta) = rownames(Omega) = rownames(Pi) = rownames(Delta) = rownames(st_B) = rownames(st_count) = rownames(Eta) = rownames(loc)
  colnames(sc_B)=colnames(sc_count)=colnames(Zeta)=colnames(B) =colnames(st_B)=colnames(st_count)=colnames(Eta)=paste0("Gene", 1:P)
  
  
  return(list(sc_count=sc_count, sc_meta = sc_meta, v = v, B = B, gamma = gamma,
              st_count = st_count, loc = loc, st_meta = st_meta, s = s,
              Omega = Omega, Delta = Delta, A = A, Pi = Pi, phi = phi, st_B = st_B, 
              Eta = Eta, sparsity = sparsity, zero = zero, psi = psi, Zeta = Zeta, sc_B=sc_B))
  
  
}


## Set the parameters for generating scRNA-seq data ############################
K = 4 # number of cell types
M = 500 # number of cells
P_k = c(0.4, 0.1, 0.2, 0.3)
N_k = M*P_k
psi_mean = 10 # means of dispersion parameter
P = 1000 # number of genes
P_gamma = 20 # number of informative genes
sigma_mean = 0.3 # mean of variance of random error 
beta_0 = 2 # baseline of log normalized expression levels
qq = 2 # common ratio or difference
e_k = c(log(1), log(3), log(6), log(18))
n_replicate = 30 # number of replicates

## Set the parameters for generating spatial data
ct.select = paste0("cellType", 1:K)
phi_mean = 10 # means of dispersion parameter
zero = 0
seed = 123

## Load location information
load("data/olfactory_bulb_patterns_new.Rdata")
st_meta = data.frame("x" = pattern_loc[, "x"], "y" = pattern_loc[, "y"], "Domain" = patterns[, "pattern_ii"]+1)
table(st_meta$Domain)
loc <- st_meta[, c("x", "y")]

## generate simulated data
for (zero in c(0.05, 0.1,0.3)) {
  for (sparsity in c(0.1, 0.3,0.5)) {
    for (replicate in 1:30) {
      seed = seed + replicate
      
      data_name = paste0("mob_pattern_zero_", zero*100, "_sparsity_", sparsity*100, "_replicate_", replicate)
      print(data_name)
      
      pseudo_data <- simulate.data(M=M, P=P, P_gamma=P_gamma, K=K, N_k=N_k, zero=zero, sparsity=sparsity, seed = seed)
        
      ## Retrieve reference data
      sc_count = pseudo_data$sc_count
      sc_meta = pseudo_data$sc_meta
      v = pseudo_data$v
      B = pseudo_data$B
      gamma = pseudo_data$gamma
      psi = pseudo_data$psi
      Zeta = pseudo_data$Zeta
      sc_B = pseudo_data$sc_B
      
      ## Retrieve ST data
      st_count = pseudo_data$st_count
      
      ## Check genes order should be matched
      if (!identical(colnames(sc_count), colnames(st_count))){
        stop("The gene names of scRNA-seq data and ST data should be matched each other!")
      }
      
      loc = pseudo_data$loc
      st_meta = pseudo_data$st_meta
      s = pseudo_data$s
      Omega = pseudo_data$Omega
      Delta = pseudo_data$Delta
      A = pseudo_data$A
      Pi = pseudo_data$Pi
      phi = pseudo_data$phi
      st_B = pseudo_data$st_B
      Eta = pseudo_data$Eta
      sparsity = pseudo_data$sparsity
      zero = pseudo_data$zero
      
      names(v) = rownames(sc_count)
      names(s) = rownames(st_count)
      
      cellType <- sort(unique(sc_meta$cellType))
      K <- length(cellType)
      
      domainType <- sort(unique(st_meta$Domain))
      D <- length(domainType)
      
      # file_name = paste0("data/simulated_data/", data_name, '.RData')
      # 
      # save(sc_count = sc_count, sc_meta = sc_meta, v = v, B = B, cellType = cellType, K = K, domainType = domainType, D = D,
      #      gamma = gamma, st_count = st_count, loc = loc, st_meta = st_meta, s = s, 
      #      Omega = Omega, Delta = Delta, A = A, Pi = Pi, phi = phi, st_B = st_B, 
      #      Eta = Eta, sparsity = sparsity, zero = zero,
      #      file = file_name)
    }
  }
}
