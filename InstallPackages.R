install.packages(c('facetscales','stringr','ggpubr','igraph','tidyverse','npreg','gridExtra','ggplot2','np','DescTools','IndepTest','mgcv','broom','furrr','RBGL','pcalg','SID','Matrix'),repos="https://mirrors.dotsrc.org/cran/")

install.packages('devtools','glmnet','mboost')
library(devtools)
install_url('https://cran.r-project.org/src/contrib/Archive/CAM/CAM_1.0.tar.gz')

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("RBGL")


devtools::install_github("zeehio/facetscales")
