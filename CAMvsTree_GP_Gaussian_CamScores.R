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
source("CAMmodified.R")

Gen_SHD_SID <- function(p,N,type){

  p <- p
  samplesize <- N
  tree_gen_type <- type
  
#Generate True Graph / Adjacency Matrix
Adjacency <- Generate_tree_adjacencymatrix(type = tree_gen_type, p= p)

Adj.Dat <- Generate_tree_adjacency_dat(Adjacency= Adjacency,p=p)

g.true.graphNEL<-as(Adjacency, "graphNEL")
g.true.igraph<-graph_from_adjacency_matrix(Adjacency)

#Generate Data


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

Data <-   Generate_GP_Data("Gaussian",noisevariances,samplesize,Adj.Dat,bandwidth=1)

#DatawoNoise <- Data[,(p+1):(2*p)]
#colnames(DatawoNoise) <- c(paste0("V",seq(1,p,1)))

Data <- Data[,1:p] 



#CAM
est.CAM <- CAM(X=Data %>% as.matrix(), scoreName = "SEMGAM", parsScore = list(numBasisFcts = 10), numCores = 1, 
               maxNumParents = 1, output = FALSE, 
               variableSel = FALSE, variableSelMethod = selGamBoost, 
               variableSelMethodPars = list(atLeastThatMuchSelected = 0.02, 
                                            atMostThatManyNeighbors = 10), pruning = FALSE, pruneMethod = selGam, 
               pruneMethodPars = list(cutOffPVal = 0.001, numBasisFcts = 10), intervData = FALSE, 
               intervMat = NA)

est.CAM.Adj <- est.CAM$Adj
g.est.CAM.graphNEL<-as(est.CAM.Adj, "graphNEL")
g.est.CAM.igraph<-graph_from_adjacency_matrix(est.CAM.Adj)


CAMScore <- ifelse(est.CAM$ScoreMat$scoreMat>-Inf, est.CAM$ScoreMat$scoreMat, 0)
EdgeScoresAsGraphNEL<-as(CAMScore, "graphNEL")
Optimal<-edmondsOptimumBranching(EdgeScoresAsGraphNEL)
AdjMatrix <- matrix(rep(0,ncol(Data)*ncol(Data)),ncol=ncol(Data))

for(i in 1:ncol(Optimal$edgeList)){
  from <- Optimal$edgeList["from",i] %>% as.integer
  to <- Optimal$edgeList["to",i] %>% as.integer
  AdjMatrix[from,to] <- 1
}

g.est.Tree.graphNEL<-as(AdjMatrix, "graphNEL")
g.est.Tree.igraph<-graph_from_adjacency_matrix(AdjMatrix)

CamScore <- {CAMScore*as(est.CAM.Adj, "matrix")} %>% rowSums() %>% sum()
TreeScore <- {CAMScore*AdjMatrix} %>% rowSums %>%  sum()

shd_tree <- shd(g.true.graphNEL,g.est.Tree.graphNEL)
shd_CAM <- shd(g.true.graphNEL,g.est.CAM.graphNEL)

sid_tree <- structIntervDist(g.true.graphNEL,g.est.Tree.graphNEL)
sid_cam <- structIntervDist(g.true.graphNEL,g.est.CAM.graphNEL)

Res <- data.frame(shd_tree = shd_tree, 
                  shd_cam = shd_CAM, 
                  sid_tree = sid_tree$sid,
                  sid_cam = sid_cam$sid,
                  sid_lower_tree = sid_tree$sidLowerBound,
                  sid_lower_cam = sid_cam$sidLowerBound,
                  sid_upper_tree = sid_tree$sidUpperBound,
                  sid_upper_cam = sid_cam$sidUpperBound,
                  CamScore = CamScore,
                  TreeScore = TreeScore
                  ) %>% 
  as_tibble(.name_repair = 'unique')
return(Res)
}

timestamp <- as.character(Sys.time()) %>% {str_replace_all(.,"[: -]","")} %>% str_sub(.,3,-3)

seed = 1
type <- c(1,2) #tree_gen_type
p <- c(16,32)
N <- c(50,100,200,500)
repetitions <- 200

sims <- expand_grid(type=type,p=p,N=N,rep=seq(1,repetitions,1))
rows <- sample(nrow(sims))
sims <- sims[rows,]


plan(multisession, workers = 60)


SummaryDat <- expand_grid(sims) %>% 
  mutate(res = future_pmap(.l=list(p=p,N=N,type=type),
                           .f=Gen_SHD_SID,
                           .progress = TRUE,
                           .options = furrr::furrr_options(seed = seed)))
future:::ClusterRegistry("stop")

saveRDS(SummaryDat,file=paste0("Data/Data_CAMvsTree_GP_Gaussian_CamScores_",timestamp,".RDS"))
