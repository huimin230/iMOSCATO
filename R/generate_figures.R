################################################################################
## 
## Title        : Robust Bayesian Integrative Modeling of Single Cell and Spatially Resolved Transcriptomics Data
## Project      : iMOSCATO
## Research goal: Bayesian cell-type deconvolution with spatial domain detection
## Authors      : Huimin Li
## Contact      : huimin.li@utdallas.edu
## Code         : Generate all figures in the paper
## Modified     : 2025-04-01
## 
################################################################################

rm(list = ls())

## Set directory
setwd("/Users/Huimin Li/PhD/Project/iMOSCATO/Paper/gitHub")
source("R/imoscato.R")


################################################################################
##
## Part 1. Simulation study
##
################################################################################

## load results
load("result/simulation_study_result.RData")

## figure 2A
p <- plot.domain(loc, size = 3, domain=Domain, colors=st_color, title = "MOB pattern")
print(p)


## figure 2B
p <- c()
sparsity_range <- c(0.1,0.3,0.5)

for (i in 1:length(prop_result)) {
  sparsity = sparsity_range[i]
  p[[i]] <- plot.pie(proportion=prop_result[[i]],
                       spatial_location=loc,
                       colors = sc_color,
                       radius = 0.5,
                       nrow=4) +
    ggtitle(paste0("Cell-ype sparsity = ", sparsity*100, "%")) +
    theme(plot.title = element_text(hjust = 0.5, size = 12)) +
    theme(legend.position = "right")
  
}

ggarrange(plotlist=p,common.legend = TRUE, legend = "right",nrow = 1)


## figure 2C
methods <- c("iMOSCATO","iMOSCATO (NB)", "CARD", "RCTD")

plot.boxplot <- function(data, ylab = "", colors = color_setting[1:length(methods)]){
  p = ggplot(data, aes(x=method, y = value, fill = method)) +
    geom_boxplot(outlier.shape = NA) +
    theme_bw() +
    facet_grid(zero~sparsity) +
    theme(legend.position="bottom") +
    labs(fill = "Method", x = "", y = ylab) +
    theme(plot.caption = element_text(hjust=0.5)) +
    guides(fill = guide_legend(nrow = 1)) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x=element_blank()) +
    scale_fill_manual(values=colors)
  
  return(p)
}

metrics <- "RMSE"
data <- metrics_result[metrics_result$metrics == metrics, ]
p <- plot.boxplot(data, ylab = metrics)
print(p)

metrics <- "PCC"
data <- metrics_result[metrics_result$metrics == metrics, ]
p <- plot.boxplot(data, ylab = metrics)
print(p)


## figure 2D
metrics <- "ARI"
data <- metrics_result[metrics_result$metrics == metrics, ]
p <- plot.boxplot(data, ylab = metrics, colors = c("#E31A1C", "#FF6EB4", "#A6CEE3", "#FDBF6F"))
print(p)


## figure S2A
metrics <- "ARI_dominant"
data <- metrics_result[metrics_result$metrics == metrics, ]
p <- plot.boxplot(data, ylab = "ARI")
print(p)


## figure S2B
data <- metrics_result[metrics_result$metrics == "AUC" & metrics_result$method == "iMOSCATO", ]
data$zero <- factor(data$zero, levels = unique(data$zero), labels = c("Low", "Medium", "High"))

p <- ggplot(data, aes(x=zero, y = value, fill = zero)) +
  geom_boxplot(outlier.shape = NA) +
  theme_bw() +
  facet_wrap(~sparsity) +
  theme(legend.position="bottom") +
  labs(fill = "Zero-inflation level", x = "", y = "AUC") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x=element_blank()) +
  scale_fill_manual(values=c("#33FF00", "#FFD92F", "#FF6EB4"))
print(p)


################################################################################
##
## Part 2. Real data analysis
##
################################################################################

datasets <- c("mouse_olfactory_bulb", "PDAC-B")
ii <- 2
data_name <- datasets[ii]

## load results
load(paste0("result/", data_name, "_result.RData"))

## figure 3A/4A manual annotation
p <- plot.domain(loc, size = 3, domain=Domain, colors=st_color, title = "Manual Annoation")
print(p)

## figure 3A/4A spatial domain for each method
methods <- c("iMOSCATO", "iMOSCATO (NB)", "BayesSpace", "Louvain")
n.methods <- length(methods)
domain_fig <- c()

for (i in 1:n.methods) {
  temp <- domain_result[domain_result$method == methods[i], ]
  temp <- temp[!is.na(temp$domain), ]
  ARI <- ari(temp$Domain, temp$domain)
  NMI <- aricode::NMI(temp[!is.na(temp$Domain), "Domain"], temp[!is.na(temp$Domain), "domain"])
  
  title = paste0(methods[i],"\n (ARI = ",round(ARI,digits = 3), ")")
  
  if (ii == 1){
    if (i %in% 1:2){
      tt <- c(1,2,4,3)
    }else if (i %in% 3:4){
      tt <- c(1,2,4,3)
    }
  }else {
    if (i %in% 1:3){
      tt <- c(2,3,1)
    }else if (i %in% 3:4){
      tt <- c(1,3,2)
    }
  }
  
  
  domain_fig[[i]] <- plot.domain(temp[, c("x", "y")],size=2.5,domain=temp$domain,
                                 label = "Domain", colors=st_color[tt], ncol = 4,title = title)+
    theme(plot.title = element_text(hjust = 0.5)) 
}

ggarrange(plotlist=domain_fig, ncol = 2,nrow=2, common.legend = TRUE, legend="bottom")


## figure 3B/4B and 3C/4C
methods <- c("iMOSCATO", "iMOSCATO (NB)", "CARD", "RCTD")
n.methods <- length(methods)
scatterpie <- c()
dominant_fig <- c()
cor_fig <- c()
type_result <- c()

for (i in 1:n.methods) {
  proportion <- prop_result[prop_result$method == methods[i], cellType]
  names(proportion) <- cellType
  spatial_location <- prop_result[prop_result$method == methods[i], c("x", "y")]
  rownames(proportion) = rownames(spatial_location) = paste0(spatial_location$x,"x",spatial_location$y)
  
  scatterpie[[i]] <- plot.pie(proportion = proportion,
                              spatial_location = spatial_location,
                              colors = sc_color[colSums(proportion>0)>0],
                              radius = 0.5,
                              nrow = 1) +
    ggtitle(methods[i]) +
    theme(plot.title = element_text(hjust = 0.5, size = 14)) 
  
  
  temp <- prop_result[!is.na(prop_result$Domain), c("Domain", "dominant", "x", "y", "method")]
  temp <- temp[temp$method == methods[i], ]
  ARI <- ari(temp$Domain, temp$dominant)
  NMI <- aricode::NMI(temp$Domain, temp$dominant)
  
  title = paste0(methods[i],"\n (ARI = ",round(ARI,digits = 3), ")")
  
  dominant_fig[[i]] <- plot.domain(temp[, c("x", "y")],size=2.5,domain=temp$dominant,
                                 label = "Domain", colors=sc_color[cellType%in%temp$dominant], ncol = 4,title = title)+
    theme(plot.title = element_text(hjust = 0.5)) 
  
  temp <- plot.prop(
    proportion = proportion,
    spatial_location = spatial_location,
    ct.visualize = cellType[-1],                 ### selected cell types to visualize
    colors = c("lightblue","lightyellow","red"), ### if not provide, we will use the default colors
    NumCols = length(cellType[-1]),
    pointSize = 1)
  
  tt <- temp$data
  tt$Method <- methods[i]
  type_result <- rbind(type_result, tt)
  
}

ggarrange(plotlist=scatterpie, ncol = 2,nrow=2, common.legend = TRUE, legend="bottom")
ggarrange(plotlist=dominant_fig, ncol = 2,nrow=2, common.legend = TRUE, legend="bottom")


## figure 3D/4D and figure 3SE/4SE
type_result$Method <- factor(type_result$Method, levels = methods)
ggplot(type_result, aes(x, y)) + 
  geom_point(aes(colour = value),size =0.5)  +
  scale_color_gradientn(colours = c("lightblue","lightyellow","red")) + 
  scale_x_discrete(expand = c(0, 1)) + scale_y_discrete(expand = c(0,1))+ 
  facet_grid(Cell_Type~Method)+
  coord_fixed()+
  theme(panel.background = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
        axis.title =element_blank())+
  labs(color = "Proportion")


## figure S3A/S4A
corrplot(cor(pcc_result), 
         method = "circle",
         type = "lower", 
         tl.col = "black", 
         tl.srt = 45, 
         tl.offset = 0.5,
         addCoef.col = "white",
         mar = c(0,0,0,0),
         diag = FALSE)


## figure S3B/S4B
plot.trace.chain(trace_result)


## figure S3C/S4C
temp <- sparsity_result$percent

p <- ggplot(temp, aes(x=method, y=value, fill=method)) +
  geom_boxplot(width = 0.5) +
  labs(x="", fill = "Method", y = "Percent of zero probability") +
  theme_bw() +
  theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size=14),
        axis.title.y = element_text(size=14),
        legend.title = element_text(size=14),
        legend.text = element_text(size=14),
        axis.text.y = element_text(size=12)) +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5, size = 14)) 
print(p)

## figure S3D/S4D
colors <- c("#6E98FF", "#7FC97F", "#E7298A", "#FFD92F", "#BEAED4", "#ff7d0b")
temp <- sparsity_result$number
p <- ggplot(data=temp, aes(x=x, y=Count, fill=x)) +
  geom_bar(stat="identity",width = 0.6) +
  xlab("Number of cell type") +
  ylab("Number of spot") +
  theme_classic() +
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        legend.position = "none") +
  scale_fill_manual(values=colors[1:length(unique(temp$Count))]) 
print(p)

