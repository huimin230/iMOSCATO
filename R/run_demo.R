################################################################################
## 
## Title        : Robust Bayesian Integrative Modeling of Single Cell and Spatially Resolved Transcriptomics Data
## Project      : iMOSCATO
## Research goal: Bayesian cell-type deconvolution with spatial domain detection
## Authors      : Huimin Li
## Contact      : huimin.li@utdallas.edu
## Code         : Tutorial to run iMOSCATO
## Modified     : 2025-04-01
## 
################################################################################


rm(list=ls())

## set your own directory
setwd("/Users/Huimin Li/PhD/Project/iMOSCATO/Paper/gitHub")
source("R/imoscato.R")

## load demo data
load("data/demo.RData")


## Create iMOSCATO object ######################################################
iMOSCATO.object <- create.iMOSCATO(
  sc_count = sc_count, 
  sc_meta = sc_meta, 
  st_count = st_count, 
  loc = loc,
  cutoff_sample = 100, 
  cutoff_feature = 0.1,
  norm_method = "tss", 
  platform = "ST")


## Run iMOSCATO ################################################################
iMOSCATO.object <- run.iMOSCATO(
  iMOSCATO.object = iMOSCATO.object, 
  n.domain = D, 
  iter = 2000,
  burn = 1000)


## Retrieve results ############################################################
## Estimated cell type proportions
prop_result <- iMOSCATO.object$result$prop_result
head(prop_result)

colors = c("#6E98FF", "#7FC97F", "#E7298A", "#FFD92F")
colnames(prop_result) <- 1:4
p = CARD::CARD.visualize.pie(proportion = prop_result, 
                   spatial_location = iMOSCATO.object$object@loc,
                   colors = colors) +
  theme(legend.title=element_text(size = 8),
        legend.text=element_text(size = 8),
        legend.box.spacing = unit(0, "pt")) 
print(p)

png("figure/imoscato_prop.png", width = 650, height = 600, res = 300)
grid.arrange(p, nrow = 1)
dev.off()


## Estimated spatial domains
domain_result = iMOSCATO.object$result$domain_result
head(domain_result)

p <- plot.domain(domain_result[,c("x","y")], size = 2, domain = domain_result$domain, colors = c("red", "steelblue3"))
print(p)

png("figure/imoscato_domain.png", width = 1000, height = 600, res = 300)
grid.arrange(p, nrow = 1)
dev.off()
