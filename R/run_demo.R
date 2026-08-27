################################################################################
## 
## Title    : A Bayesian Integrative Model of Single-cell and Spatially Resolved Transcriptomics Data
## Project  : iMOSCATO
## Authors  : Huimin Li
## Contact  : huimin.li01@utrgv.edu
## Modified : 2025-08-20
## 
################################################################################


rm(list=ls())

## set your own directory
setwd("/Users/huiminli/Desktop/iMOSCATO")
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
p = plot.pie(proportion = prop_result, 
                   spatial_location = iMOSCATO.object$object@loc,
                   colors = colors,
             nrow = 1) +
  theme(legend.title=element_text(size = 8),
        legend.text=element_text(size = 8),
        legend.box.spacing = unit(0, "pt")) 
print(p)

png("figure/imoscato_prop.png", width = 700, height = 600, res = 300)
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
