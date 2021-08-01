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

Gen_IC_BER <- function(p,N,type){
  
samplesize <- N
tree_gen_type <- type
  
#Generate True Graph / Adjacency Matrix
Adjacency <- Generate_tree_adjacencymatrix(type = tree_gen_type, p= p)
#Generate True Graph / Adjacency Matrix
Adj.Dat <- Generate_tree_adjacency_dat(Adjacency= Adjacency,p=p)

g.true.graphNEL<-as(Adjacency, "graphNEL")
g.true.igraph<-graph_from_adjacency_matrix(Adjacency)

#Generate Data



#Generating Data

noisevariances <- c(runif(1,min=1,max=2),runif(p-1,min=1/25,max=2/25))

CausalFunctions <- Generate_CausalFunctions(Adj.Dat=Adj.Dat,bandwidth=1)
Data <- Generate_Data_From_Causal_Functions("Gaussian",noisevariances,samplesize,Adj.Dat,bandwidth=1,p,CausalFunctions)

Eds <- Edmonds(Data,noise="G",workers=1)


#Calculate True Gaussian Score of graph  
##### CHANGE THIS SUCH THAT WE JUST SUM NORMALIZED EDGE WEIGHTS INSTEAD OF THIS COMPLICATED
# THAT IS EQUIVALENT.
CalcScore <- function(Eds,EdgeList){

#Calc AdjMatrix from EdgeList    
AdjMatrix <- matrix(rep(0,p*p),ncol=p)
  for(i in 1:ncol(EdgeList)){
    from <- EdgeList["from",i] %>% as.integer
    to <- EdgeList["to",i] %>% as.integer
    AdjMatrix[from,to] <- 1
  }
#Find Rootnode
rootnode <- which(colSums(AdjMatrix)==0)

#Calc Adjacency Dataframe
Adj.Dat <- data.frame(parent=0,node = rootnode) 
for(i in 1:(p-1)){
  for(j in (i+1):p){
    if(AdjMatrix[i,j]!=0)
    {
      Adj.Dat <- rbind(Adj.Dat,data.frame(parent=i,node = j)  )
    }
  }}
Adj.Dat <- Adj.Dat %>% arrange(node)

#Score is TargetWeight of Root + non-normalized EdgeWeight on edges
SumOfNonNormalizedEdgeWeights <- left_join(Adj.Dat,Eds$EdgeWeights,
                                           by=c("parent"="from","node"="to"))$EdgeWeight %>% 
  sum(log(.),na.rm =TRUE)
RootWeight <- Eds$TargetWeights %>% 
  filter(node==rootnode) %>% 
  select(TargetWeight) %>% 
  pull %>% 
  log()

Score <- SumOfNonNormalizedEdgeWeights+RootWeight
return(Score)
}


BestScore <- CalcScore(Eds,Eds$Optimal$edgeList)
#Find Second best scoring graph:
Scores <- rep(-Inf,p-1)
for( i in 1:p-1){
  from <- Eds$Optimal$edgeList["from",i] %>% as.numeric()
  to <- Eds$Optimal$edgeList["to",i] %>% as.numeric()
  NewEdgeWeightMatrix <- Eds$EdgeWeightMatrix
  NewEdgeWeightMatrix[from,to] <- 0
  rownames(NewEdgeWeightMatrix) <- sub(".", "", rownames(NewEdgeWeightMatrix))
  colnames(NewEdgeWeightMatrix) <- sub(".", "", colnames(NewEdgeWeightMatrix))
  NewEdgeScoresAsGraphNEL<-as(NewEdgeWeightMatrix, "graphNEL")
  Optimal<-edmondsOptimumBranching(NewEdgeScoresAsGraphNEL)
  Scores[i] <-CalcScore(Eds,Optimal$edgeList)
}
SecondBestScore <- min(Scores)

IdentifiabilityConstant <- 1/2*(SecondBestScore-BestScore) #We May divide to 1/2

CalcEdgeReversal <- function(from,to){
  R <- gam(Data[,paste0("V",as.character(from))] %>% as.matrix~s(Data[,paste0("V",as.character(to))] %>% as.matrix,bs='tp'))$residuals
  HR <- KLentropy(R,k=10)$Estimate
  Hto <-  KLentropy(Data[,paste0("V",as.character(to))] %>% as.matrix,k=10)$Estimate
  HRto <- KLentropy(cbind(R,Data[,paste0("V",as.character(to))] %>% as.matrix),k=10)$Estimate
  ER <- HR + Hto - HRto
  return(ER)
}
CalcEdgeReversal <- Vectorize(CalcEdgeReversal)

MinimumEdgeReversal <- Adj.Dat %>% 
  filter(parent != 0) %>% 
  mutate(ER= CalcEdgeReversal(from=parent,to=node)) %>% 
  select(ER) %>% 
  pull %>% 
  min()

# V1 <- Data[,"V1"] %>% as.matrix()
# V2 <- Data[,"V2"] %>% as.matrix()
# V3 <- Data[,"V3"] %>% as.matrix() 
# V4 <- Data[,"V4"] %>% as.matrix()
# 
# H1 <-   KLentropy(V1,k=10)$Estimate
# H12 <-  KLentropy(cbind(V1,V2),k=10)$Estimate
# H13 <-  KLentropy(cbind(V1,V3),k=10)$Estimate
# H123 <- KLentropy(cbind(V1,V2,V3),k=10)$Estimate
# 
# H2 <-   KLentropy(V2,k=10)$Estimate
# H23 <-  KLentropy(cbind(V2,V3),k=10)$Estimate
# H24 <-  KLentropy(cbind(V2,V4),k=10)$Estimate
# H234 <- KLentropy(cbind(V2,V3,V4),k=10)$Estimate
# 
# CI1 <- H12 + H13 - H123 - H1
# CI2 <- H23 + H24 - H234 - H2
return(data.frame(IC=IdentifiabilityConstant, MER= MinimumEdgeReversal))
}


timestamp <- as.character(Sys.time()) %>% {str_replace_all(.,"[: -]","")} %>% str_sub(.,3,-3)
 
seed = 1
type <- c(1) #tree_gen_type
p <- c(8,16)
N <- c(200000)
repetitions <- 100



sims <- expand_grid(type=type,p=p,N=N,repetitions=seq(1,repetitions,1))
rows <- sample(nrow(sims))
sims <- sims[rows,]

#plan(multisession, workers = 4)
plan(multicore, workers = 60)


SummaryDat <- expand_grid(sims) %>% 
  mutate(res = future_pmap(.l=list(p=p,N=N,type=type),
                           .f=Gen_IC_BER,
                           .progress = TRUE,
                           .options = furrr::furrr_options(seed = seed)))
future:::ClusterRegistry("stop")

saveRDS(SummaryDat,file=paste0("Data/Data_IdentifiabilityConstantMultivariate_",timestamp,".RDS"))
