library(igraph)
library(tidyverse)
library(npreg)
library(gridExtra)
library(ggplot2)
library(DescTools)
library(IndepTest)
library(mgcv)
library(broom)
library(furrr)
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

Gen_SHD_SID <- function(p,N,a){
  p <- p
  samplesize <- N
  alpha <- a
  
#Generate True Graph / Adjacency Matrix
Adjacency <- matrix(0,ncol=p,nrow=p)

for(i in 1:(p-1)){
  for(j in (i+1):p){
    if(sum(Adjacency[,j])==0){    
      Adjacency[i,j] <- sample(c(0,1),size=1,prob=c(0.9,0.1))}
    else {
      Adjacency[i,j] <- 0}
    if(j==i+1 & sum(Adjacency[,j])==0){
      Adjacency[i,j] <- 1
    }
  }
}

Adj.Dat <- data.frame(parent=0,node = 1) 

for(i in 1:(p-1)){
  for(j in (i+1):p){
    if(Adjacency[i,j]!=0)
    {
      Adj.Dat <- rbind(Adj.Dat,data.frame(parent=i,node = j)  )
    }
  }}
Adj.Dat <- Adj.Dat %>% arrange(node)


g.true.graphNEL<-as(Adjacency, "graphNEL")
g.true.igraph<-graph_from_adjacency_matrix(Adjacency)

#Generate Data

#Function to generate gaussian process causal functions
Generate_Causal_Functions <- function(p,samplesize,Adj.Dat,bandwidth){
  l <- bandwidth #bandwidth?
  Data <- matrix(0,ncol=p,nrow=samplesize) 
  Data[,1] <- mgcv:::rmvn(1, rep(0,samplesize), diag(rep(1,samplesize)))
  
  for(nod in 2:p){
    d = abs(outer(Data[,Adj.Dat[nod,1]],Data[,Adj.Dat[nod,1]],"-")) # compute distance matrix, d_{ij} = |x_i - x_j|
    Sigma_SE = exp(-d^2/(2*l^2)) # squared exponential kernel
    Data[,nod] <- mgcv:::rmvn(1, rep(0,samplesize), Sigma_SE)
  }
  colnames(Data) <- paste0("V",seq(1,p,1),sep="")
  Data <- Data %>% as_tibble()
  
  RegressionDirections <- Adj.Dat %>% filter(parent != 0) %>% rename(from=parent,to=node) %>% as_tibble()
  
  Res <- RegressionDirections %>% 
    mutate(f = list(purrr:::pmap(.l=list(from,to),.f=function(from,to){ 
      form <- formula(paste0("V",to,"~","s(V",from,",bs='tp')"))
      return(gam(form, data = Data))
    }))) 
  
  return(Res)
}

#Function to generate noise innovations:

Generate_noise <- function(samplesize,noisedist,noisevariances,p){
  Noise <- matrix(0,nrow=samplesize,ncol=p) 
  if(noisedist=="Gaussian"){
    for(i in 1:p){
      Noise[,i] <- rnorm(samplesize,mean=0,sd=sqrt(noisevariances[i]))
    }
  } else if(noisedist == "t"){
    for(i in 1:p){
      Noise[,i] <- rt(samplesize, df=4, ncp=0) #If Noise is dist=t then noisevariances ignored var ~1
    }
  }
  else if(is.numeric(noisedist)){
    for(i in 1:p){
      temp <- rnorm(samplesize,mean=0,sd=sqrt(noisevariances[i]))
      Noise[,i] <- sign(temp)*abs(temp)^noisedist
    }
  }
  Noise
}


# Combining GP causal function and noise innovation generation:
Generate_GP_Data <- function(noisedist,noisevariances,samplesize,Adj.Dat,bandwidth){
  
  Noise <- Generate_noise(samplesize,noisedist,noisevariances,p)
  
  l <- bandwidth #bandwidth?
  Data <- matrix(0,ncol=p,nrow=samplesize) 
  DatawoNoise <- matrix(0,ncol=p,nrow=samplesize) 
  Data[,1] <- Noise[,1]
  DatawoNoise[,1] <- Data[,1]
  for(nod in 2:p){
    d = abs(outer(Data[,Adj.Dat[nod,1]],Data[,Adj.Dat[nod,1]],"-")) # compute distance matrix, d_{ij} = |x_i - x_j|
    Sigma_SE = exp(-d^2/(2*l^2)) # squared exponential kernel
    DatawoNoise[,nod] <- mgcv:::rmvn(1, rep(0,samplesize), Sigma_SE)
    Data[,nod] <- DatawoNoise[,nod] + Noise[,nod]
  }
  
 #Data <- cbind(Data,DatawoNoise) %>% as_tibble(.name_repair = NULL)
  colnames(Data) <- paste0("V",seq(1,p,1),sep="")
  Data <- Data %>% as_tibble()
  return(Data)
}



#Generating Data
noisevariances <- c(runif(1,min=1,max=2),runif(p-1,min=1/25,max=2/25))
Data <-   Generate_GP_Data(noisedist=alpha,noisevariances,samplesize,Adj.Dat,bandwidth=1)
Data <- Data[,1:p] 


#CAT
ptm <- proc.time()
est.CAT <- CAT(Data,noise="G",workers=1)
time.ectg <- proc.time() - ptm

g.est.TreeG.graphNEL<-as(est.CAT$AdjMatrix, "graphNEL")
g.est.TreeG.igraph<-graph_from_adjacency_matrix(est.CAT$AdjMatrix)

ptm <- proc.time()
est.CAT <- CAT(Data,noise="NA",workers=1)
time.ecte <- proc.time() - ptm

g.est.TreeNG.graphNEL<-as(est.CAT$AdjMatrix, "graphNEL")
g.est.TreeNG.igraph<-graph_from_adjacency_matrix(est.CAT$AdjMatrix)

#CAM
ptm <- proc.time()
est.CAM <- CAM(X=Data %>% as.matrix(), scoreName = "SEMGAM", parsScore = list(numBasisFcts = 10), numCores = 1, 
               maxNumParents = 1, output = FALSE, 
               variableSel = FALSE, variableSelMethod = selGamBoost, 
               variableSelMethodPars = list(atLeastThatMuchSelected = 0.02, 
                                            atMostThatManyNeighbors = 10), pruning = FALSE, pruneMethod = selGam, 
               pruneMethodPars = list(cutOffPVal = 0.001, numBasisFcts = 10), intervData = FALSE, 
               intervMat = NA)
time.cam <- proc.time() - ptm

est.CAM.Adj <- est.CAM$Adj
g.est.CAM.graphNEL<-as(est.CAM.Adj, "graphNEL")
g.est.CAM.igraph<-graph_from_adjacency_matrix(est.CAM.Adj)


#SHD
shd_treeG <- shd(g.true.graphNEL,g.est.TreeG.graphNEL)
shd_treeNG <- shd(g.true.graphNEL,g.est.TreeNG.graphNEL)
shd_CAM <- shd(g.true.graphNEL,g.est.CAM.graphNEL)

#SID
sid_treeG <- structIntervDist(g.true.graphNEL,g.est.TreeG.graphNEL)
sid_treeNG <- structIntervDist(g.true.graphNEL,g.est.TreeNG.graphNEL)
sid_cam <- structIntervDist(g.true.graphNEL,g.est.CAM.graphNEL)


Res <- data.frame(shd_treeG = shd_treeG, 
                  shd_treeNG = shd_treeNG, 
                  shd_cam = shd_CAM, 
                  sid_treeG = sid_treeG$sid,
                  sid_treeNG = sid_treeNG$sid,
                  sid_cam = sid_cam$sid,
                  sid_lower_treeG = sid_treeG$sidLowerBound,
                  sid_lower_treeNG = sid_treeNG$sidLowerBound,
                  sid_lower_cam = sid_cam$sidLowerBound,
                  sid_upper_treeG = sid_treeG$sidUpperBound,
                  sid_upper_treeNG = sid_treeNG$sidUpperBound,
                  sid_upper_cam = sid_cam$sidUpperBound,
                  time.ectg = time.ectg["elapsed"],
                  time.ecte = time.ecte["elapsed"],
                  time.cam = time.cam["elapsed"]) %>% 
  as_tibble(.name_repair = 'unique')
return(Res)
}

timestamp <- as.character(Sys.time()) %>% {str_replace_all(.,"[: -]","")} %>% str_sub(.,3,-3)

seed = 1
set.seed(1)

p <- c(32)
N <- c(50,500)
a <- c(seq(0.1,2,0.1),seq(2.5,4,0.5))
repetitions <- 500

sims <- expand_grid(p=p,N=N,a=a,rep=seq(1,repetitions,1))
rows <- sample(nrow(sims))
sims <- sims[rows,]

plan(multisession, workers = 60)


SummaryDat <-  sims %>% 
  mutate(res = future_pmap(.l=list(p=p,N=N,a=a),
                           .f=Gen_SHD_SID,
                           .progress = TRUE,
                           .options = furrr::furrr_options(seed = seed)))
future:::ClusterRegistry("stop")

saveRDS(SummaryDat,file=paste0("Data/Data_CAMvsTree_GP_NonGaussian_",timestamp,".RDS"))
