library(igraph)
library(tidyverse)
library(npreg)
library(gridExtra)
library(ggplot2)
library(ggpubr)
library(facetscales)
library(grid)
#library(np)
library(DescTools)
library(IndepTest)
library(mgcv)
library(broom)
library(furrr)
library(stringr)
library(RBGL)
library(pcalg) # SHD
library(SID)  # SID
#install.packages("glmnet") #needed for CAM  
#install.packages("mboost") #needed for CAM
#install.packages("CAM_1.0.tar.gz", repos = NULL, type="source")
library(CAM)
library(Matrix)

source("Functions.R")


#Generating Data
set.seed(1000)

p <- list()
for(i in 2:5){
  l <- 1
  x <- seq(-2,2,0.01)
  d = abs(outer(x,x,"-")) # compute distance matrix, d_{ij} = |x_i - x_j|
  Sigma_SE = exp(-d^2/(2*l^2)) # squared exponential kernel
  res <- mgcv:::rmvn(1, rep(0,length(x)), Sigma_SE)
  dat <- data.frame(x=x,y=res)
  p[[i-1]]<- ggplot(data=dat)+
    geom_line(aes(x=x,y=y))+
    theme(axis.title.x = element_blank())+
    theme(axis.title.y = element_blank())
}

Plot <- ggarrange(p[[1]],p[[2]],p[[3]],p[[4]], nrow=1, common.legend = TRUE, legend="right")
Plot <- annotate_figure(Plot, left = textGrob("f(x)", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                bottom = textGrob("x", gp = gpar(cex = 1.3)))


ggsave(
  "Plots/GPsamplepaths.jpeg",
  plot = Plot,
  device = "jpeg",
  path = NULL,
  scale = 1,
  width = 210,
  height = 75,
  units = c("mm"),
  dpi = 500,
  limitsize = TRUE
)
                