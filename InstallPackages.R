install.packages(c('hash','facetscales','stringr','ggpubr','igraph','tidyverse','npreg','gridExtra','ggplot2','np','DescTools','IndepTest','mgcv','broom','furrr','RBGL','pcalg','SID','Matrix','Rcpp','readxl','gnlearn','bnlearn','xtable','kableExtra','latex2exp','glmnet','mboost'),lib="/home/mj/R/x86_64-pc-linux-gnu-library/4.1",repos="https://mirrors.dotsrc.org/cran/")

install.packages('devtools',lib="/home/mj/R/x86_64-pc-linux-gnu-library/4.1")
library(devtools)
install_url('https://cran.r-project.org/src/contrib/Archive/CAM/CAM_1.0.tar.gz')

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("RBGL")

install.packages("https://www.bnlearn.com/releases/bnlearn_latest.tar.gz", repos = NULL, type = "source")

if (!requireNamespace('devtools', quietly=TRUE))
  install.packages('devtools')
devtools::install_github('rlebron-bioinfo/gnlearn')

devtools::install_github("zeehio/facetscales")
