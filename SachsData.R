library(igraph)
library(tidyverse)
library(ggplot2)
library(IndepTest)
library(mgcv)
library(broom)
library(RBGL)
library(pcalg) # SHD
library(SID)  # SID
library(CAM)
library(Matrix)
library(Rcpp)
library("readxl")
library(gnlearn)
library(bnlearn)
library(ggpubr)
library(expm)
library(xtable)
library(kableExtra)
library(latex2exp)

source("Functions.R")
source("CAT.R")


addsmall <-function(x){return(x+rnorm(length(x),mean=0,sd=0.000001))}
std <-function(x){return((x-mean(x))/sd(x))}


data <- bind_rows(read_xls(path="SachsData/1. cd3cd28.xls")) %>%  
  mutate_all(.f=addsmall) %>%   
  mutate_all(.funs=std)

#CAT 
G <- CAT(data,noise="G",workers=1,numbasisfnct = 25)
NG <- CAT(data,noise="NA",workers=1,numbasisfnct = 25)
#CAT.Prune
G001 <- CAT(data,noise="G",workers=1,numbasisfnct = 25, pvalCutoff = 0.001)
NG001 <- CAT(data,noise="NA",workers=1,numbasisfnct = 25, pvalCutoff = 0.001)
#CAM
p <- 4
cam.est <- CAM(X=data %>% as.matrix(), scoreName = "SEMGAM", parsScore = list(numBasisFcts = 25), numCores = 1, 
    maxNumParents = p, output = FALSE, 
    variableSel = TRUE, variableSelMethod = selGamBoost, 
    variableSelMethodPars = list(atLeastThatMuchSelected = 0.02, 
                                 atMostThatManyNeighbors = 10), 
    pruning = TRUE, pruneMethod = selGam, 
    pruneMethodPars = list(cutOffPVal = 0.001, numBasisFcts = 10), intervData = FALSE, 
    intervMat = NA)


#GES
score <- new("GaussL0penObsScore", data)
ges.fit <- pcalg::ges(score, phase=c('forward','backward'), iterate=TRUE)

#notears 
notearsAdjMat <- gnlearn::notears(  data,  lambda1 = 0.1,  loss.type = c("l2"),  max.iter = 100,  h.tol = 1e-8,  rho.max = 1e+16,
  w.threshold = 0.1,  m = NULL,  to = c("adjacency"),  seed = 1) %>%  apply(.,2,FUN=function(x){ifelse(x>0,1,0)})

#MMHC
#mmhc bge
mmhcbge.fit <- bnlearn::mmhc(data, whitelist = NULL, blacklist = NULL, restrict.args = list(),
     maximize.args = list(score="bge"), debug = FALSE)
#mmhc bic
mmhcbic.fit <- bnlearn::mmhc(data, whitelist = NULL, blacklist = NULL, restrict.args = list(),
                    maximize.args = list(score="bic-g"), debug = FALSE)

#Create estimated adjacency matrices
GAdjMat <- apply(G$EdgeWeightMatrix*as.matrix(G$AdjMatrix),2,FUN=function(x){ifelse(x>0,1,0)})
NGAdjMat <- apply(NG$EdgeWeightMatrix*as.matrix(NG$AdjMatrix),2,FUN=function(x){ifelse(x>0,1,0)})
G001AdjMat <- apply(G001$EdgeWeightMatrix*as.matrix(G001$AdjMatrix),2,FUN=function(x){ifelse(x>0,1,0)})
NG001AdjMat <- apply(NG001$EdgeWeightMatrix*as.matrix(NG001$AdjMatrix),2,FUN=function(x){ifelse(x>0,1,0)})
camAdjMat <- NGAdjMat*0 + cam.est$Adj
gesAdjMat <- as(as(ges.fit$repr,"graphNEL"),"matrix")
mmhcbgeAdjMat <- GAdjMat*matrix(rep(0,11*11),ncol=11)
for (n in 1:nrow(mmhcbge.fit$arcs)){
  mmhcbgeAdjMat[mmhcbge.fit$arcs[n,1],mmhcbge.fit$arc[n,2]] <- 1
}
mmhcbicAdjMat <- GAdjMat*matrix(rep(0,11*11),ncol=11)
for (n in 1:nrow(mmhcbic.fit$arcs)){
  mmhcbicAdjMat[mmhcbic.fit$arcs[n,1],mmhcbic.fit$arc[n,2]] <- 1
}

#Create true adjacency matrix
truth <- data.frame(from=c("praf","pmek",  "plcg","PIP3","PIP3","PKA","PKA",     "PKA",   "PKA", "PKA", "PKA", "PKC","PKC","PKC", "PKC",  "PKC",   "p44/42","PIP3","plcg","PIP2"),
                    to=c(  "pmek","p44/42","PIP2","plcg","PIP2","P38","pakts473","p44/42","pmek","praf","pjnk","P38","PKA","praf","pjnk", "pmek", "pakts473","pakts473","PKC","PKC")) %>% 
  mutate_all(.f=as.character)
#Creating adjacency matrix for true graph
TrueAdjMat <- GAdjMat*matrix(rep(0,11*11),ncol=11)
for (n in 1:nrow(truth)){
  TrueAdjMat[truth[n,1],truth[n,2]] <- 1
}


#emptyAdjMat
emptyAdjMat <- GAdjMat*matrix(rep(0,11*11),ncol=11)


# g.CATG.igraph<-graph_from_adjacency_matrix(GAdjMat)
# g.CATNG.igraph<-graph_from_adjacency_matrix(NGAdjMat)
# g.CATG001.igraph<-graph_from_adjacency_matrix(G001AdjMat)
# g.CATNG001.igraph<-graph_from_adjacency_matrix(NG001AdjMat)
# g.CAM.igraph<-graph_from_adjacency_matrix(camAdjMat)
# g.true.igraph<-graph_from_adjacency_matrix(TrueAdjMat)
# g.GES.igraph <- graph_from_adjacency_matrix(gesAdjMat)
# g.notears.igraph <- graph_from_adjacency_matrix(notearsAdjMat)
# g.mmhcbge.igraph <- graph_from_adjacency_matrix(mmhcbgeAdjMat)
# g.mmhcbic.igraph <- graph_from_adjacency_matrix(mmhcbicAdjMat)
# 
# par(mfrow=c(3,3))
# 
# plot(g.true.igraph,main="Truth")
# plot(g.CATG.igraph,main="CAT.G")
# plot(g.CATNG.igraph,main="CAT.E")
# plot(g.CATG001.igraph,main="CAT.G")
# plot(g.CATNG001.igraph,main="CAT.E")
# plot(g.CAM.igraph,main="CAM.G")
# plot(g.GES.igraph,main="GES BIC")
# plot(g.notears.igraph,main="NoTears")
# plot(g.mmhcbge.igraph,main="MMHC BGe")
# plot(g.mmhcbic.igraph,main="MMHC BIC")

TrueAdjData <- Generate_tree_adjacency_dat(TrueAdjMat,11)

Data <- tibble(Method = c("CAT", "CAT","CAT", "CAT", "CAM", "GES","NoTears","MMHC", "MMHC","EmptyGraph"),
               Score = c("Entropy","Gaussian","Entropy","Gaussian","Gaussian","BIC","---","BGe","BIC","---"),
               Prune =c("No","No","Yes","Yes","Yes","---","---","---","---","---"),
  AdjecencyMatrix = list(NGAdjMat, GAdjMat,NG001AdjMat, G001AdjMat,camAdjMat, gesAdjMat,notearsAdjMat,mmhcbgeAdjMat,mmhcbicAdjMat,emptyAdjMat),
               AdjacencyData= purrr::map(.x= AdjecencyMatrix,.f=function(x){Generate_tree_adjacency_dat(x,11)}),
               TrueAdjDat=list(TrueAdjData),
               TrueAdjMat = list(TrueAdjMat))

#Random single rooted graph generation function
GenerateRandomSingleRootedDag <- function(p){
  adj <- Generate_tree_adjacencymatrix(type = 1, p= p) %>% as.matrix()
  #Add extra edges
  adj[which(adj==0)] <- sample(c(0,1),size=length(which(adj==0)),replace=TRUE,prob=c(0.95,0.05))
  adj[upper.tri(adj)==FALSE] <- 0
  return(adj)
}

RandomGraphData <- tibble(rep = seq(1,100,1),
       Method = "RandomGraph",
       Score = "---",
       Prune = "---",
       p = 11,
       AdjecencyMatrix = map(.x=p,.f=GenerateRandomSingleRootedDag),
       AdjacencyData= purrr::map(.x= AdjecencyMatrix,.f=function(x){Generate_tree_adjacency_dat(x,11)}),
       TrueAdjDat=list(TrueAdjData),
       TrueAdjMat = list(TrueAdjMat))



AnalyzeDataFixedp <- function(Data,p){
  
  Res <- Data %>% 
    mutate(
      SHD = pmap_dbl(.l=list(AdjecencyMatrix,TrueAdjMat),.f=function(x,y){
        EstGraph = as(x %>% as.matrix(), "graphNEL")
        TrueGraph = as(y %>% as.matrix(),"graphNEL")
        pcalg::shd(EstGraph,TrueGraph)
      }),
      SHD.CPDAG = pmap_dbl(.l=list(AdjecencyMatrix,TrueAdjMat),.f=function(x,y){
        EstGraph = as(x %>% as.matrix(), "graphNEL")
        TrueGraph = as(y %>% as.matrix(),"graphNEL")
        pcalg::shd(pcalg::dag2cpdag(EstGraph),pcalg::dag2cpdag(TrueGraph))
      }),
      SID = pmap_dbl(.l=list(AdjecencyMatrix,TrueAdjMat),.f=function(x,y){
        EstGraph = as(x, "graphNEL")
        TrueGraph = as(y, "graphNEL")
        structIntervDist(TrueGraph,EstGraph)$sid
      }),
      CorrectPredEdges = pmap_dbl(.l=list(AdjacencyData,TrueAdjDat),.f=function(x,y){
        CorrectEdges <- inner_join(x,y,by=c("parent","node"))
        n <- nrow(CorrectEdges)
        n
      }),
      TotalPredEdges = pmap_dbl(.l=list(AdjacencyData),.f=function(x){
        TotalEdges <- x
        n <- nrow(TotalEdges)
        n
      }),
      TotalTrueEdges = pmap_dbl(.l=list(TrueAdjDat),.f=function(x){
        TotalEdges <- x
        n <- nrow(TotalEdges)
        n
      }),
      Ed.CorrectPredOverTotalPred = CorrectPredEdges/TotalPredEdges,
      Ed.CorrectPredOverTotalTrue = CorrectPredEdges/TotalTrueEdges,
      Anc.MeanNodeWiseCorrectPredOverTotalPred = pmap_dbl(.l=list(AdjecencyMatrix,TrueAdjMat,p),.f=function(x,y,p){
        Est.Ancestors <- Find_Ancestors(x,p)
        True.Ancestors <- Find_Ancestors(y,p)
        PercentCorrect <- NA
        for(i in 1:p){
          x <- True.Ancestors[i,2] %>% pull %>% unlist() %>% as.numeric()
          y <- Est.Ancestors[i,2] %>% pull %>% unlist() %>% as.numeric()
          if(length(y)==0){
            PercentCorrect[i] <- 1
          } else{
            PercentCorrect[i] <- length(intersect(x,y))/length(y)
          }
        }
        mean(PercentCorrect)
      }),
      Anc.CorrectPredOverTotalPred = pmap_dbl(.l=list(AdjecencyMatrix,TrueAdjMat,p),.f=function(x,y,p){
        Est.Ancestors <- Find_Ancestors(x,p)
        True.Ancestors <- Find_Ancestors(y,p)
        Correct <- NA
        Total <- NA
        for(i in 1:p){
          x <- True.Ancestors[i,2] %>% pull %>% unlist() %>% as.numeric()
          y <- Est.Ancestors[i,2] %>% pull %>% unlist() %>% as.numeric()
          Correct[i] <- length(intersect(x,y))
          Total[i] <- length(y)
          
        }
        sum(Correct)/sum(Total)
      }),
      Anc.CorrectPredOverTrue = pmap_dbl(.l=list(AdjecencyMatrix,TrueAdjMat,p),.f=function(x,y,p){
        Est.Ancestors <- Find_Ancestors(x,p)
        True.Ancestors <- Find_Ancestors(y,p)
        Correct <- NA
        Total <- NA
        for(i in 1:p){
          x <- True.Ancestors[i,2] %>% pull %>% unlist() %>% as.numeric()
          y <- Est.Ancestors[i,2] %>% pull %>% unlist() %>% as.numeric()
          Correct[i] <- length(intersect(x,y))
          Total[i] <- length(x)
          
        }
        sum(Correct)/sum(Total)
      })
    )
  
}  


resRandomData <- AnalyzeDataFixedp(RandomGraphData,11) %>% 
  group_by(Method,Score) %>% 
  summarize(SHD=mean(SHD),SHD.CPDAG=mean(SHD.CPDAG),SID=mean(SID),Ed.CorrectPredOverTotalPred=mean(Ed.CorrectPredOverTotalPred),Ed.CorrectPredOverTotalTrue = mean(Ed.CorrectPredOverTotalTrue))


res <- AnalyzeDataFixedp(Data,11)

res %>% 
  select(Method,Prune,Score,SHD,SHD.CPDAG,SID,Ed.CorrectPredOverTotalPred,Ed.CorrectPredOverTotalTrue,CorrectPred= CorrectPredEdges,TotalPred= TotalPredEdges,TotalTrue= TotalTrueEdges ) %>% 
  arrange(SID,SHD)


options(knitr.kable.NA = '---')

res %>% 
  select(Method,Prune,Score,SHD,SHD.CPDAG,SID,Ed.CorrectPredOverTotalPred,Ed.CorrectPredOverTotalTrue) %>% 
  bind_rows(resRandomData) %>% 
  arrange(SHD,SID) %>% 
  kable(format = 'latex', booktabs = TRUE,digits=3)  

#### Testing

vars <- c(truth %>% select(from) %>%  pull %>%  as.character(), truth %>%  select(to) %>%  pull %>% as.character()) %>% unique()

Hyps <- expand_grid(from = vars, to = vars) %>% 
  filter(from != to) %>% 
  mutate(hyp = 1)
Hyps <- bind_rows(Hyps,Hyps %>% mutate(hyp=-1))

translation <- data.frame(original = as.character(colnames(data)),new=seq(1,ncol(data),1))
Hyps <- Hyps %>% mutate(from=as.character(from),to=as.character(to)) %>% 
  rowwise() %>% 
  mutate(from = translation[which(translation[,1]==from),2],
         to = translation[which(translation[,1]==to),2]) %>% 
  ungroup()



source("Functions.R")
TestRes <- HypothesisTest(dat = data,hyp=Hyps,alpha = 0.05,numbasisfnct = 25)

saveRDS(TestRes,file="Data/TestRes.RDS")
TestRes <- readRDS(file="Data/TestRes.RDS")

TestRes <- TestRes %>% 
  rowwise() %>% 
  mutate(from = translation[which(translation[,2]==from),1],
         to = translation[which(translation[,2]==to),1]) %>% 
  ungroup() 



Hypotheses <- TestRes %>% 
  left_join(truth %>% mutate(IsEdgePresent = TRUE),by = c("from"="from","to"="to")) %>% 
  mutate(IsEdgePresent = ifelse(is.na(IsEdgePresent),FALSE,IsEdgePresent),
         dist = ifelse(IsEdgePresent == TRUE, 1, 0))

AdjacencyDist <- matrix(rep(0,11*11),ncol=11)
k <- 1
while(k <= 11-1){
  AdjacencyDist <- AdjacencyDist + (TrueAdjMat%^%k)*k
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

stat_total <- Hypotheses %>% 
  mutate(dist = ifelse(dist >= 0 , ifelse(dist== 0,dist,1),-1)) %>% 
  group_by(hyp,IsEdgePresent) %>% 
  summarize(#TestExact = mean(Test_Exact_Bonferroni),
    TestAsymp=mean(Test_Asym_Bonferroni)) %>% 
  gather(c(TestAsymp),key=Type,val=Test) %>% 
  mutate(Stat = ifelse( (IsEdgePresent == "TRUE" & hyp == 1) | (IsEdgePresent == FALSE & hyp == -1 ),"Level","Power")) %>% 
  mutate(col = paste(Stat,IsEdgePresent,sep="_")) %>% 
  ungroup %>% 
  select(-IsEdgePresent,-Stat,-hyp) %>% 
  arrange(col) %>% 
  pivot_wider(names_from = col, values_from = Test)

stat_total_n <- Hypotheses %>% 
  mutate(dist = ifelse(dist >= 0 , ifelse(dist== 0,dist,1),-1)) %>% 
  group_by(hyp,IsEdgePresent) %>% 
  summarize(n=n()) %>% 
  gather(c(n),key=Type,val=Test) %>% 
  mutate(Stat = ifelse( (IsEdgePresent == "TRUE" & hyp == 1) | (IsEdgePresent == FALSE & hyp == -1 ),"Level","Power")) %>% 
  mutate(col = paste(Stat,IsEdgePresent,sep="_")) %>% 
  ungroup %>% 
  select(-IsEdgePresent,-Stat,-hyp) %>% 
  pivot_wider(names_from = col, values_from = Test)

stat_total <- bind_rows(stat_total,stat_total_n)

Power_bydist <- Hypotheses %>% 
  mutate(dist = ifelse(dist >= 0 , ifelse(dist== 0,dist,1),-1)) %>% 
  group_by(hyp,IsEdgePresent,dist) %>% 
  summarize(#TestExact = mean(Test_Exact_Bonferroni),
    TestAsymp=mean(Test_Asym_Bonferroni)) %>% 
  gather(c(TestAsymp),key=Type,val=Test) %>% 
  mutate(Stat = ifelse( (IsEdgePresent == "TRUE" & hyp == 1) | (IsEdgePresent == FALSE & hyp == -1 ),"Level","Power")) %>% 
  filter(Stat == "Power") %>% 
  mutate(col = paste(Stat,IsEdgePresent,dist,sep="_")) %>% 
  ungroup %>% 
  select(-IsEdgePresent,-Stat,-hyp,-dist) %>% 
  arrange(col) %>% 
  pivot_wider(names_from = col, values_from = Test)

Power_bydist_n <- Hypotheses %>% 
  mutate(dist = ifelse(dist >= 0 , ifelse(dist== 0,dist,1),-1)) %>% 
  group_by(hyp,IsEdgePresent,dist) %>% 
  summarize(n = n()) %>% 
  gather(c(n),key=Type,val=Test) %>% 
  mutate(Stat = ifelse( (IsEdgePresent == "TRUE" & hyp == 1) | (IsEdgePresent == FALSE & hyp == -1 ),"Level","Power")) %>% 
  filter(Stat == "Power") %>% 
  mutate(col = paste(Stat,IsEdgePresent,dist,sep="_")) %>% 
  ungroup %>% 
  select(-IsEdgePresent,-Stat,-hyp,-dist) %>% 
  arrange(col) %>% 
  pivot_wider(names_from = col, values_from = Test)

Power_bydist <- bind_rows(Power_bydist,Power_bydist_n)

Table <- left_join(stat_total,Power_bydist) %>% 
  select(Type,'Power_FALSE_-1','Power_FALSE_1','Power_FALSE_0','Power_FALSE','Power_TRUE','Level_TRUE','Level_FALSE')


#PrintTableToLatex:
options(knitr.kable.NA = '---')

Table %>%  
  mutate_at(.vars=c('Power_FALSE_-1','Power_FALSE_1','Power_FALSE_0','Power_FALSE','Power_TRUE'),.funs=function(x){floor(x* 100) / 100}) %>% 
  mutate_at(.vars=c('Level_TRUE','Level_FALSE'),.funs=function(x){ceiling(x* 100) / 100}) %>%   
  kable(format = 'latex', booktabs = TRUE) %>% 
  add_header_above(header = c("0" = 2, "Power" = 5,"Level"=2))
