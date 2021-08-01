library(igraph)
library(tidyverse)
library(npreg)
library(gridExtra)
library(ggplot2)
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

setwd("/home/lnd974/trees/Simulations")

source("Functions.R")


Gen_Data <- function(lambda,alpha,N){

  Xt <- rnorm(N,mean=0,sd=1)
  X <- sign(Xt)*abs(Xt)^alpha
  Y <-  (1-lambda)*(X^3)+lambda*X + rnorm(N,mean=0,sd=1)
  Data <- cbind(X,Y)
  return(Data)
}

Gen_IdentifiabilityConstant <- function(lambda,alpha,N){
  
    Data <- Gen_Data(lambda=lambda,alpha=alpha,N=N)
    
    X <- Data[,"X"]
    Y <- Data[,"Y"]
    R <- gam(X~s(Y,bs='tp',k=10))$residuals
    
    HRY <- KLentropy(cbind(R,Y),k=10)$Estimate
    HR <- KLentropy(R,k=10)$Estimate
    HY <- KLentropy(Y,k=10)$Estimate
    
    MI <- HR + HY - HRY
    pval <- MINTauto(R,Y,kmax=10,B1=1000,B2=1000)[1]
    
    return(data.frame(MI=MI,pval=pval))
}


timestamp <- as.character(Sys.time()) %>% {str_replace_all(.,"[: -]","")} %>% str_sub(.,3,-3)

seed = 1
alpha <- seq(0.25,1.75,0.1)
lambda <- seq(0,1,0.05)
N <- c(50000)


sims <- expand_grid(alpha,lambda,N)
rows <- sample(nrow(sims))
sims <- sims[rows,]

plan(multisession, workers = 60)


SummaryDat <- sims %>% 
  mutate(res = future_pmap(.l=list(lambda=lambda,alpha=alpha,N=N),
                           .f=Gen_IdentifiabilityConstant,
                           .progress = TRUE,
                           .options = furrr::furrr_options(seed = seed)))
future:::ClusterRegistry("stop")

saveRDS(SummaryDat,file=paste0("Data/Data_IdentifiabilityConstant_",timestamp,".RDS"))