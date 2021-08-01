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


Gen_SHD_SID <- function(p,N,type){
  p <- p
  samplesize <- N
  tree_gen_type <- type
  
#Generate Tree Graph / Adjacency Matrix
Adjacency <- Generate_tree_adjacencymatrix(type = tree_gen_type, p= p) %>% as.matrix()

#Add extra edges
Adjacency[which(Adjacency==0)] <- sample(c(0,1),size=length(which(Adjacency==0)),replace=TRUE,prob=c(0.95,0.05))
Adjacency[upper.tri(Adjacency)==FALSE] <- 0

Adj.Dat.True <- Generate_tree_adjacency_dat(Adjacency= Adjacency,p=p)

Ancestors <- Find_Ancestors(Adjacency,p)


g.true.graphNEL<-as(Adjacency, "graphNEL")
g.true.igraph<-graph_from_adjacency_matrix(Adjacency)

#Generate Data
noisevariances <- c(runif(1,min=1,max=2),runif(p-1,min=1/25,max=5/25))
CausalFunctions <- Generate_CausalFunctions(Adj.Dat=Adj.Dat.True,bandwidth=1)
Data <- Generate_Data_From_Causal_Functions("Gaussian",noisevariances,samplesize,Adj.Dat.True,bandwidth=1,p,CausalFunctions)


start_time = Sys.time()
est.CAT <- CAT(Data,noise="G",workers=1)
end_time = Sys.time()
CAT_Time <- difftime(end_time,start_time,units="secs") %>% as.numeric()

start_time = Sys.time()
est.CAM.VarSel <- CAM(X=Data %>% as.matrix(), scoreName = "SEMGAM", parsScore = list(numBasisFcts = 10), numCores = 1, 
                     maxNumParents = p, output = FALSE, 
                     variableSel = TRUE, variableSelMethod = selGamBoost, 
                     variableSelMethodPars = list(atLeastThatMuchSelected = 0.02, 
                                                  atMostThatManyNeighbors = 10), 
                     pruning = FALSE, pruneMethod = selGam, 
                     pruneMethodPars = list(cutOffPVal = 0.001, numBasisFcts = 10), intervData = FALSE, 
                     intervMat = NA)
end_time = Sys.time()
CAM.VarSel_Time <- difftime(end_time,start_time,units="secs") %>% as.numeric()

start_time = Sys.time()
est.CAM.PruneVarSel <- CAM(X=Data %>% as.matrix(), scoreName = "SEMGAM", parsScore = list(numBasisFcts = 10), numCores = 1, 
                           maxNumParents = p, output = FALSE, 
                           variableSel = TRUE, variableSelMethod = selGamBoost, 
                           variableSelMethodPars = list(atLeastThatMuchSelected = 0.02, 
                                                        atMostThatManyNeighbors = 10), 
                           pruning = TRUE, pruneMethod = selGam, 
                           pruneMethodPars = list(cutOffPVal = 0.001, numBasisFcts = 10), intervData = FALSE, 
                           intervMat = NA)
end_time = Sys.time()
CAM.PruneVarSel_Time <- difftime(end_time,start_time,units="secs") %>% as.numeric()


Adj.Data <- tibble(Method = c("True","CAT.G","CAM.VarSel","CAM.PruneVarSel"),
       EstimationTime = c(0,CAT_Time,CAM.VarSel_Time,CAM.PruneVarSel_Time),
       AdjecencyMatrix = list(Adjacency %>%  as.matrix(),
                  est.CAT$AdjMatrix %>% as.matrix(),
                  est.CAM.VarSel$Adj %>% as.matrix(),
                  est.CAM.PruneVarSel$Adj %>% as.matrix())) %>% 
  mutate(AdjacencyData = pmap(.l=list(AdjecencyMatrix),.f=function(AdjecencyMatrix){
    Generate_tree_adjacency_dat(AdjecencyMatrix,p)
  }))


return(Adj.Data)
}

timestamp <- as.character(Sys.time()) %>% {str_replace_all(.,"[: -]","")} %>% str_sub(.,3,-3)

seed = 1
type <- c(1) #tree_gen_type
p <- c(16,32,64)
N <- c(50,250,500)
repetitions <- 100


sims <- expand_grid(type=type,p=p,N=N,rep=seq(1,repetitions,1))
rows <- sample(nrow(sims))
sims <- sims[rows,]


plan(multicore, workers = 60)


SummaryDat <- expand_grid(sims) %>% 
  mutate(res = future_pmap(.l=list(p=p,N=N,type=type),
                           .f=Gen_SHD_SID,
                           .progress = TRUE,
                           .options = furrr::furrr_options(seed = seed)))
future:::ClusterRegistry("stop")

saveRDS(SummaryDat,file=paste0("Data/Data_CAMvsTree_GP_Gaussian_SingleRootedDags_",timestamp,".RDS"))
