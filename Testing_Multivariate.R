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
library(stringr)
library(RBGL)
library(pcalg) # SHD
library(SID)  # SID
#library(CAM)
library(Matrix)
library(expm)

setwd("/home/lnd974/trees/Simulations")

source("Functions.R")


Gen_Test <- function(p,N,TreeType,level){
  p <- p
  samplesize <- N
  tree_gen_type <- TreeType
  
  #Generate Tree Graph / Adjacency Matrix
  Adjacency <- Generate_tree_adjacencymatrix(type = tree_gen_type, p= p) %>% as.matrix()
  
  Adj.Dat.True <- Generate_tree_adjacency_dat(Adjacency= Adjacency,p=p)
  
  Hypotheses <- expand_grid(from=seq(1,p,1),to=seq(1,p,1)) %>% 
    filter(from != to) %>% 
    left_join(Adj.Dat.True %>% mutate(IsEdgePresent = TRUE),by = c("from"="parent","to"="node")) %>% 
    mutate(IsEdgePresent = ifelse(is.na(IsEdgePresent),FALSE,IsEdgePresent),
           dist = ifelse(IsEdgePresent == TRUE, 1, 0))
  
  AdjacencyDist <- matrix(rep(0,p*p),ncol=p)
  k <- 1
  while(k <= p-1){
    AdjacencyDist <- AdjacencyDist + (Adjacency%^%k)*k
    k <- k +1
  }

  for (i in 1:nrow(Hypotheses)){
    if (Hypotheses[i,"IsEdgePresent"] == TRUE){
      next
    } 
    from = Hypotheses[i,"from"] %>% pull
    to = Hypotheses[i,"to"] %>% pull
    if (AdjacencyDist[from,to] > 0){
      Hypotheses[i,"dist"] <- AdjacencyDist[from,to]
    } else if (AdjacencyDist[to,from] > 0){
      Hypotheses[i,"dist"] <- -AdjacencyDist[to,from]
    } 
  }

  Hypothesis <- bind_rows(Hypotheses %>% select(from,to) %>%  mutate(Hyp = 1),Hypotheses %>% select(from,to) %>%  mutate(Hyp = -1))
  #Generate Data
  noisevariances <- c(runif(1,min=1,max=2),runif(p-1,min=1/(25),max=5/(25)))
  CausalFunctions <- Generate_CausalFunctions(Adj.Dat=Adj.Dat.True,bandwidth=1)
  dat <- Generate_Data_From_Causal_Functions("Gaussian",noisevariances,samplesize,Adj.Dat.True,bandwidth=1,p,CausalFunctions)
  
  Test = HypothesisTest(dat = dat, hyp = Hypothesis,alpha= level, numbasisfnct = NULL)
  Results = Hypotheses %>% left_join(Test,by=c("from"="from","to"="to"))
  return(Results)
}

timestamp <- as.character(Sys.time()) %>% {str_replace_all(.,"[: -]","")} %>% str_sub(.,3,-3)

set.seed(1)
seed = 1
TreeType <- c(2) #tree_gen_type
p <- c(2,4,6,8,16)
N <- c(500,1000,5000,10000,20000)
level <- c(0.05)
repetitions <- c(400)


sims <- expand_grid(TreeType=TreeType,p=p,N=N,rep=seq(1,repetitions,1),level = level)
rows <- sample(nrow(sims))
sims <- sims[rows,]


plan(multisession, workers = 60)

SummaryDat <- expand_grid(sims) %>% 
  mutate(Test = future_pmap(.l=list(p=p,N=N,TreeType=TreeType,level=level),
                           .f=Gen_Test,
                           .progress = TRUE,
                           .options = furrr::furrr_options(seed = seed)))
future:::ClusterRegistry("stop")

message("Parallization code has been executed without error")

DataNameWOtimestamp <- paste0("Data_Testing_Multivariate_")

saveRDS(SummaryDat,file=paste0("Data/",DataNameWOtimestamp,timestamp,".RDS"))

message(paste0("Data-save has been executed without error: ",DataNameWOtimestamp,timestamp,".RDS"))

file.copy(from="Testing_Multivariate.R", to=paste0("Data/Log/",DataNameWOtimestamp,timestamp,".R"), overwrite = TRUE, copy.mode = TRUE, copy.date = FALSE)

message("Log-save has been executed without error")