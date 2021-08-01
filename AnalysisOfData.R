library(igraph)
library(tidyverse)
library(npreg)
library(gridExtra)
library(ggplot2)
library(ggpubr)
library(facetscales)
library(DescTools)
library(IndepTest)
library(mgcv)
library(broom)
library(furrr)
library(stringr)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("RBGL")
library(RBGL)
library(pcalg) # SHD
library(SID)  # SID
#install.packages("glmnet") #needed for CAM  
#install.packages("mboost") #needed for CAM
#install.packages("CAM_1.0.tar.gz", repos = NULL, type="source")
library(CAM)
library(Matrix)
library(viridis)

source("Functions.R")


############################
# Identifiability Constant #
############################


Data <- readRDS(file="Data/Data_IdentifiabilityConstant_2105051615.RDS")  %>%  
  unnest(cols="res") %>% 
  filter(alpha >= 0.25) %>% 
  filter(alpha <=1.75)

Data %>% select(alpha) %>% unique()

MI <- ggplot() +
  geom_raster(data=Data,aes(x=lambda,y=alpha,fill=MI))+
  geom_tile(data=Data %>% filter(pval>=0.05),aes(x=lambda,y=alpha,fill=MI,color="red"),size=0.6,linetype=1)+
  ylab(expression(paste("Non-Gaussianity: ",alpha)))+
  xlab(expression(paste("Linearity: ",lambda)))+
  scale_fill_viridis(option = "D")+
  labs(fill = "Identifiability\nGap")+
  scale_colour_manual(name = expression("p-value">="0.05"),
                      values =c('red'='red'), labels = c(''),
                      guide = guide_legend(override.aes = list(fill="white") ))+
  theme(panel.background = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "cm"),
        
        panel.spacing = unit(0, "cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=1))+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0))


ggsave(
  "Plots/IdentifiabilityConstant_Heatmap.pdf",
  plot = MI,
  device = cairo_pdf,
  path = NULL,
  scale = 1,
  width = 210,
  height = 80,
  units = c("mm"),
  dpi = 500,
  limitsize = TRUE
)

#########################################
# Identifiability Constant Multivariate #
#########################################

#Read Data
Data <- readRDS(file="Data/Data_IdentifiabilityConstantMultivariate_2105120958.RDS")  %>%  unnest(cols="res")

plevel <- Data %>% select(p) %>% unique() %>% arrange(p) %>%  pull(p) %>% paste()

Data <- Data %>%  mutate(p= factor(p,levels = plevel))

ICMulti <- Data %>% mutate(Difference = IC-MER) %>% 
  ggplot(data=.,aes(x=p,y=Difference,group=p)) +
  geom_hline(yintercept=0,color="red",linetype=2)+
  #geom_violin()+
  geom_boxplot(width=0.1,outlier.shape = NA)+
  geom_jitter(alpha=0.4,color="deepskyblue2",shape=16,position = position_jitter(height = 0, width = .4))+
  scale_y_continuous(limits=c(-0.11,0.25))+
  ylab("Identifiability Gap \n minus minimum edge reversal")


ggsave(
  "Plots/IdentifiabilityConstant_Multivariate.pdf",
  plot = ICMulti,
  device = cairo_pdf,
  path = NULL,
  scale = 1,
  width = 210,
  height = 75,
  units = c("mm"),
  dpi = 500,
  limitsize = TRUE
)

#Procent of models with positive difference:
Data %>% mutate(Difference = IC-MER) %>% 
  mutate(hh = ifelse(Difference>0,1,0)) %>% 
  group_by(hh) %>% 
  summarise(n=n())


#################
# GRAPH EXAMPLE #
#################
# Plotting two randomly generated Type1 and Type2 trees.

set.seed(56)

pdf(file="Plots/TreeTypes.pdf",
    width=32, height=16)
par(mfrow=c(1,2),mar=c(0,0,0,0))
p <- 100


# Type1 tree
type <- 1

Adjacency <- Generate_tree_adjacencymatrix(type = type, p= p)
Adj.Dat <- Generate_tree_adjacency_dat(Adjacency,p)

g.true.graphNEL<-as(Adjacency, "graphNEL")
g.true.igraph<-graph_from_adjacency_matrix(Adjacency)
l <- layout_with_fr(g.true.igraph)

leafnode <- ifelse(rowSums(Adjacency)>0,0,1)
sum(leafnode)
V(g.true.igraph)$color <- ifelse(leafnode==1, "chartreuse3",'darkgoldenrod3')
V(g.true.igraph)$color[1] <- 'gray29'


plot(g.true.igraph, 
     vertex.size=5, 
     vertex.label=NA,
     layout=l,
     edge.lty=1, 
     edge.arrow.size=1,
     edge.width=1,
     edge.color='black',
     arrow.color='black',
     vertex.color = V(g.true.igraph)$color) 

title("Type 1",
      adj= 0.3,
      line= -5,
      cex.main=4,
      col.main="black")

#Type 2 tree
type <- 2

Adjacency <- Generate_tree_adjacencymatrix(type = type, p= p)
Adj.Dat <- Generate_tree_adjacency_dat(Adjacency,p)

g.true.graphNEL<-as(Adjacency, "graphNEL")
g.true.igraph<-graph_from_adjacency_matrix(Adjacency)
l <- layout_with_fr(g.true.igraph)

leafnode <- ifelse(rowSums(Adjacency)>0,0,1)
sum(leafnode)
V(g.true.igraph)$color <- ifelse(leafnode==1, "chartreuse3",'darkgoldenrod3')
V(g.true.igraph)$color[1] <- 'gray29'

 plot(g.true.igraph, 
     vertex.size=5, 
     vertex.label=NA,
     layout=l,
     edge.lty=1, 
     edge.arrow.size=1,
     edge.width=1,
     edge.color='black',
     arrow.color='black',
     vertex.color = V(g.true.igraph)$color) 

 title("Type 2",
       adj= 0.3,
       line= -5,
       cex.main=4,
       col.main="black")
dev.off()

#####################
#### GP Gaussian ####
#####################

#Read Data
Data <- readRDS(file="Data/Data_CAMvsTree_GP_Gaussian_2104262319.RDS")  %>%  unnest(cols="res")


Nlevel <- Data %>% select(N) %>% unique() %>% arrange(N) %>%  pull(N) %>% paste()
typelevel <- Data %>% select(type) %>% unique() %>% arrange(type) %>%  pull(type) %>% paste()
DATA_res <- Data %>% 
  mutate(shd_diff = shd_tree-shd_cam,
         sid_diff = sid_tree-sid_cam,
         SampleSize= factor(N,levels = Nlevel),
         Type= factor(type,levels = typelevel))

#SHD FOR EACH METHOD
shd_type1 <- DATA_res %>% select(Type,p,SampleSize,rep,shd_tree,shd_cam) %>% 
  filter(Type ==1) %>% 
  rename(CAT.G=shd_tree,CAM=shd_cam) %>% 
  gather(Method,value,c(CAT.G,CAM)) %>% 
  ggplot(data=.,aes(x=SampleSize,y=value,fill=Method))+
  geom_boxplot(outlier.shape=23,
               outlier.size=1)+ 
  ylab("")+
  xlab("")+
  facet_wrap(Type ~p, scales = "free",labeller = label_both,ncol=4)+ 
  theme(legend.position="right")

shd_type2 <- DATA_res %>% select(Type,p,SampleSize,rep,shd_tree,shd_cam) %>% 
  filter(Type ==2) %>% 
  rename(CAT.G=shd_tree,CAM=shd_cam) %>% 
  gather(Method,value,c(CAT.G,CAM)) %>% 
  ggplot(data=.,aes(x=SampleSize,y=value,fill=Method))+
  geom_boxplot(outlier.shape=23,
               outlier.size=1) +
  ylab("")+
  xlab("")+
  facet_wrap(Type ~p, scales = "free",labeller = label_both,ncol=4)+ 
  theme(legend.position="right")


Plot <- ggarrange(shd_type1, shd_type2, nrow=2, common.legend = TRUE, legend="right")

Plot <- annotate_figure(Plot,
                bottom = text_grob("Sample Size", color = "black",vjust=-1),
                left = text_grob("SHD to true graph", color = "black",vjust=2, rot = 90)
)


ggsave(
  "Plots/CAMvsTree_GP_Gaussian_SHD_Type1and2.pdf",
  plot = Plot,
  device = "pdf",
  path = NULL,
  scale = 1,
  width = 210,
  height = 140,
  units = c("mm"),
  dpi = 500,
  limitsize = TRUE
)


#SID FOR EACH METHOD
sid_type1 <- DATA_res %>% select(Type,p,SampleSize,rep,sid_tree,sid_cam) %>% 
  filter(Type ==1) %>% 
  rename(CAT.G=sid_tree,CAM=sid_cam) %>% 
  gather(Method,value,c(CAT.G,CAM)) %>% 
  ggplot(data=.,aes(x=SampleSize,y=value,fill=Method))+
  geom_boxplot(outlier.shape=23,
               outlier.size=1) +
  ylab("")+
  xlab("")+
  facet_wrap(Type ~p, scales = "free",labeller = label_both,ncol=4)+ 
  theme(legend.position="right")

sid_type2 <- DATA_res %>% select(Type,p,SampleSize,rep,sid_tree,sid_cam) %>% 
  filter(Type ==2) %>% 
  rename(CAT.G=sid_tree,CAM=sid_cam) %>% 
  gather(Method,value,c(CAT.G,CAM)) %>% 
  ggplot(data=.,aes(x=SampleSize,y=value,fill=Method))+
  geom_boxplot(outlier.shape=23,
               outlier.size=1) +
  ylab("")+
  xlab("")+
  facet_wrap(Type ~p, scales = "free",labeller = label_both,ncol=4)+ 
  theme(legend.position="right")


Plot <- ggarrange(sid_type1, sid_type2, nrow=2, common.legend = TRUE, legend="right")

Plot <- annotate_figure(Plot,
                        bottom = text_grob("Sample Size", color = "black",vjust=-1),
                        left = text_grob("SID to true graph", color = "black",vjust=2, rot = 90)
)


ggsave(
  "Plots/CAMvsTree_GP_Gaussian_SID_Type1and2.pdf",
  plot = Plot,
  device = "pdf",
  path = NULL,
  scale = 1,
  width = 210,
  height = 140,
  units = c("mm"),
  dpi = 500,
  limitsize = TRUE
)


#Dataframe median

DATA_res %>% select(type,p,SampleSize,rep,shd_tree,shd_cam) %>% 
  rename(ECT.G=shd_tree,CAM=shd_cam) %>% 
  gather(Method,value,c(ECT.G,CAM)) %>% 
  group_by(type,p,SampleSize,Method) %>%
  summarise(median=median(value),mean=mean(value)) %>% 
  arrange(p,SampleSize,Method,type) %>%
  print(n=100)




### CAM SCORES CHECK ####
#Here we use the CAM scores on Edmonds algorithm to check that previous observations
# of the difference in performance is not due to possibly different conditional mean
# estimation techniques.

#Read GP Gaussian data
Data <- readRDS(file="Data/Data_CAMvsTree_GP_Gaussian_CamScores_2105042122.RDS")  %>%  unnest(cols="res")


Nlevel <- Data %>% select(N) %>% unique() %>% arrange(N) %>%  pull(N) %>% paste()
typelevel <- Data %>% select(type) %>% unique() %>% arrange(type) %>%  pull(type) %>% paste()
DATA_res <- Data %>% 
  mutate(shd_diff = shd_tree-shd_cam,
         sid_diff = sid_tree-sid_cam,
         SampleSize= factor(N,levels = Nlevel),
         Type= factor(type,levels = typelevel))

#SHD FOR EACH METHOD
shd <- DATA_res %>% select(Type,p,SampleSize,rep,shd_tree,shd_cam) %>% 
  rename(CAT.G=shd_tree,CAM=shd_cam) %>% 
  gather(Method,value,c(CAT.G,CAM)) %>% 
  ggplot(data=.,aes(x=SampleSize,y=value,fill=Method))+
  geom_boxplot(outlier.shape=23,
               outlier.size=1) +
  facet_wrap(Type ~ p, scales = "free",labeller = label_both,ncol=4)+ 
  theme(legend.position="right")+
  xlab("Sample Size")+
  ylab("SHD to true graph")


ggsave(
  "Plots/CAMvsTree_GP_Gaussian_CamScores.pdf",
  plot = shd,
  device = "pdf",
  path = NULL,
  scale = 1,
  width = 210,
  height = 75,
  units = c("mm"),
  dpi = 500,
  limitsize = TRUE
)


########################
#### GP NonGaussian ####
########################
#Read Data
Data <- readRDS(file="Data/Data_CAMvsTree_GP_NonGaussian_2104290812.RDS")  %>%  unnest(cols="res") 


alevel <- Data %>% select(a) %>%  unique() %>% pull(a) %>%  sort() %>% paste()
Nlevel <- Data %>% select(N) %>%  unique() %>% pull(N) %>%  sort() %>% paste()


DATA_res <- Data %>% 
  mutate(shd_diffG = shd_treeG-shd_cam,
         sid_diffG = sid_treeG-sid_cam,
         shd_diffNG = shd_treeNG-shd_cam,
         sid_diffNG = sid_treeNG-sid_cam,
         SampleSize= factor(N,levels = Nlevel),
         a = factor(a,levels=alevel))

#SHD FOR EACH METHOD

scales_y <- list(
  `50` = scale_y_continuous(limits = c(0, 30)),
  `500` = scale_y_continuous(limits = c(0, 15)),
  `1000` = scale_y_continuous(limits = c(0, 10))
)

shd <- DATA_res %>% select(p,a,SampleSize,rep,shd_treeG,shd_treeNG,shd_cam) %>% 
  rename(CAT.G=shd_treeG,CAT.E = shd_treeNG,CAM=shd_cam) %>% 
  gather(Method,value,c(CAT.G,CAT.E,CAM)) %>% 
  group_by(p,a,SampleSize,Method) %>% 
  summarize(mean=median(value),q20=quantile(value,probs=0.25),q80=quantile(value,probs=0.75)) %>% 
  ungroup %>% 
  ggplot(data=., aes(x=a, y=mean, colour=Method,group=Method)) +  
  geom_line() +
  geom_ribbon(aes(ymin=q20, ymax=q80,fill=Method), linetype=2, alpha=0.05)+
  facet_grid_sc(rows = vars(SampleSize), scales = list(y = scales_y),labeller = label_both)+
  ylab("SHD to true graph")+
  theme(legend.position = "right")+
  xlab(expression(alpha)) 

ggsave(
  "Plots/CAMvsTree_GP_NonGaussian.pdf",
  plot = shd,
  device = cairo_pdf,
  path = NULL,
  scale = 1,
  width = 210,
  height = 120,
  units = c("mm"),
  dpi = 500,
  limitsize = TRUE
)



#SID FOR EACH METHOD
# scales_y <- list(
#   `50` = scale_y_continuous(limits = c(0, 425)),
#   `500` = scale_y_continuous(limits = c(0, 300))
# )
# 
# sid <- DATA_res %>% select(p,a,SampleSize,rep,sid_treeG,sid_treeNG,sid_cam) %>% 
#   rename(CAT.G=sid_treeG,CAT.E = sid_treeNG,CAM=sid_cam) %>% 
#   gather(Method,value,c(CAT.G,CAT.E,CAM)) %>% 
#   group_by(p,a,SampleSize,Method) %>% 
#   summarize(median=median(value),q20=quantile(value,probs=0.25),q80=quantile(value,probs=0.75)) %>% 
#   ungroup %>% 
#   ggplot(data=., aes(x=a, y=median, colour=Method,group=Method)) + 
#   geom_point() + 
#   geom_line() +
#   geom_ribbon(aes(ymin=q20, ymax=q80,fill=Method), linetype=2, alpha=0.05)+
#   facet_grid_sc(rows = vars(SampleSize), scales = list(y = scales_y),labeller = label_both)+
#   ylab("SID to true graph")+
#   xlab(expression(alpha))
# 
# 
# Plot <- ggarrange(shd, sid, nrow=1, common.legend = TRUE, legend="right")



######################
####   TIMINGS    ####
######################

#Small sample timings SingleRootedDags:
Data <- readRDS(file="Data/Data_CAMvsTree_GP_Gaussian_SingleRootedDags_2105181606.RDS")  %>%  unnest(cols="res") 

Data %>% 
  select(type,p,N,rep,Method,EstimationTime) %>% 
  ggplot(data=.,aes(x=p,y=EstimationTime,group=interaction(N,Method,p)))+
  geom_boxplot(aes(fill=Method))

Data %>% 
  select(type,p,N,rep,Method,EstimationTime) %>% 
  group_by(type,p,N,Method) %>% 
  summarise(meantime= mean(EstimationTime))

######################
#### 3 Node setup ####
######################


Data <- readRDS(file="Data/Data_CAMvsTree_Custom3node_Gaussian_2104202046.RDS")  %>%  
  unnest(cols="res") %>%  
  select(-EdsWMG.X1,-EdsWMG.X2,-EdsWMG.X3,-EdsWMNG.X1,-EdsWMNG.X2,-EdsWMNG.X3) %>%  
  unique()



Nlevel <- c("50","100","200","500","1000","5000","10000","20000")
DATA_res <- Data %>% 
  mutate(shd_diffG = shd_treeG-shd_cam,
         sid_diffG = sid_treeG-sid_cam,
         shd_diffNG = shd_treeNG-shd_cam,
         sid_diffNG = sid_treeNG-sid_cam,
         SampleSize= factor(N,levels = Nlevel))

#SHD FOR EACH METHOD

plot <- DATA_res %>% select(SampleSize,rep,shd_treeG,shd_treeNG,shd_cam) %>% 
  group_by(SampleSize) %>%
  rename(CAT.G=shd_treeG,CAT.E = shd_treeNG,CAM=shd_cam) %>% 
  gather(Method,SHD,c("CAT.G","CAT.E","CAM")) %>% 
  group_by(SampleSize,Method) %>% 
  summarise(mean=mean(SHD),median=median(SHD)) %>% 
  ungroup %>% 
  gather(Stat.,SHD,c(mean,median)) %>% 
  ggplot(data=.,aes(x=SampleSize,y=SHD,group=interaction(Method,Stat.),color=Method,linetype=Stat.))+
  geom_line()+
  ylab("SHD to true graph")+
  xlab("Sample size")


ggsave(
  "Plots/CAMvsTree_3nodeCustom_SHD.pdf",
  plot = plot,
  device = "pdf",
  path = NULL,
  scale = 1,
  width = 210,
  height = 65,
  units = c("mm"),
  dpi = 500,
  limitsize = TRUE
)



############################
#### Single-rooted DAGS ####
############################

#Read Data
Data <- readRDS(file="Data/Data_CAMvsTree_GP_Gaussian_SingleRootedDags_2105102144.RDS")  %>%  unnest(cols="res") 


Data <- left_join(Data %>% 
                    filter(Method != "True"),
                  Data %>% 
                    filter(Method == "True") %>%  
                    rename(TrueAdjMat = AdjecencyMatrix,TrueAdjDat=AdjacencyData) %>% 
                    select(-Method,-EstimationTime),
                  by=c("type","p","N","rep"))

  



AnalyzeDataFixedp <- function(Data,k){
  
  Res <- Data %>% 
    filter(p==k) %>% 
    mutate(
      SHD = pmap_dbl(.l=list(AdjecencyMatrix,TrueAdjMat),.f=function(x,y){
        EstGraph = as(x %>% as.matrix(), "graphNEL")
        TrueGraph = as(y %>% as.matrix(),"graphNEL")
        shd(EstGraph,TrueGraph)
      }),
      SID = pmap_dbl(.l=list(AdjecencyMatrix,TrueAdjMat),.f=function(x,y){
        EstGraph = as(x, "graphNEL")
        TrueGraph = as(y, "graphNEL")
        structIntervDist(EstGraph,TrueGraph)$sid
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


Data16 <- AnalyzeDataFixedp(Data,16)
Data32 <- AnalyzeDataFixedp(Data,32)
Data64 <- AnalyzeDataFixedp(Data,64)

Data <- bind_rows(Data16,Data32,Data64)

plevel <- Data %>% select(p) %>%  unique() %>% pull(p) %>%  sort() %>% paste()
Nlevel <- Data %>% select(N) %>%  unique() %>% pull(N) %>%  sort() %>% paste()

Res <- Data %>% 
  mutate(p = factor(p,level=plevel),
         N = factor(N,level=Nlevel)) %>% 
  filter(Method != "CAM.VarSel")

Res <- Res %>% 
  mutate(Method = ifelse(Method == "CAM.PruneVarSel","CAM",Method))
#SHD

SHD <- Res %>%  
  select(p,N,rep,Method,SHD) %>% 
  ggplot(data=.,aes(x=N,y=SHD,fill=Method))+
  geom_boxplot(outlier.shape=23,
               outlier.size=1) +
  facet_wrap(~p, ncol=3,scales="free_y",labeller="label_both")+
  xlab("Sample Size")

ggsave(
  "Plots/CAMvsTree_Singlerooted_SHD.pdf",
  plot = SHD,
  device = "pdf",
  path = NULL,
  scale = 1,
  width = 210,
  height = 75,
  units = c("mm"),
  dpi = 500,
  limitsize = TRUE
)

# #SID 
# 
# Res %>%  
#   select(p,N,rep,Method,SID) %>% 
#   ggplot(data=.,aes(x=N,y=SID,color=Method))+
#   geom_boxplot()+
#   facet_wrap(~p, ncol=3,scales="free_y")

### EDGES STATISTICS ###

#Percent of Predicted edges that are correct -  Ed.CorrectPredOverTotalPred

EDGE1 <- Res %>%  
  select(p,N,rep,Method,Ed.CorrectPredOverTotalPred) %>% 
  ggplot(data=.,aes(x=N,y=Ed.CorrectPredOverTotalPred,fill=Method))+
  geom_boxplot(outlier.shape=23,
                outlier.size=1) +
  facet_wrap(~p, ncol=3,labeller="label_both")+
  xlab("")+
  ylab(expression("Edges: True positive rate"))+
  theme(plot.margin = unit(c(0,0,0,5), "mm"),
        axis.title.y= element_text(hjust = 0.5))

#Correct Predicted edges Over Total True Edges -  Ed.CorrectPredOverTotalTrue

EDGE2 <- Res %>%  
  select(p,N,rep,Method,Ed.CorrectPredOverTotalTrue) %>% 
  ggplot(data=.,aes(x=N,y=Ed.CorrectPredOverTotalTrue,fill=Method))+
  geom_boxplot(outlier.shape=23,
               outlier.size=1) +
  facet_wrap(~p, ncol=3,labeller="label_both")+
  xlab("Sample Size")+
  ylab(expression('Correctly recovered edges \nover total causal edges'))+
  theme(plot.margin = unit(c(0,0,0,5), "mm"))


Plot <- ggarrange(EDGE1, EDGE2, nrow=2, common.legend = TRUE, legend="right")

ggsave(
  "Plots/CAMvsTree_Singlerooted_Edges.pdf",
  plot = Plot,
  device = "pdf",
  path = NULL,
  scale = 1,
  width = 210,
  height = 150,
  units = c("mm"),
  dpi = 500,
  limitsize = TRUE
)


### ANCESTOR STATISTICS ###

# sum Correct Ancestors over sum of all predicted ancestors  -  Anc.CorrectPredOverTotalPred

ANC1 <- Res %>%  
  select(p,N,rep,Method,Anc.CorrectPredOverTotalPred) %>% 
  ggplot(data=.,aes(x=N,y=Anc.CorrectPredOverTotalPred,fill=Method))+
  geom_boxplot(outlier.shape=23,
               outlier.size=1) +
  facet_wrap(~p, ncol=3,labeller="label_both")+
  xlab("")+
  ylab(expression("Ancestors: True positive rate"))+
  theme(plot.margin = unit(c(0,0,0,5), "mm"))

# sum Correct Ancestors over sum of all true ancestors  -  Anc.CorrectPredOverTrue

ANC2 <- Res %>%  
  select(p,N,rep,Method,Anc.CorrectPredOverTrue) %>% 
  ggplot(data=.,aes(x=N,y=Anc.CorrectPredOverTrue,fill=Method))+
  geom_boxplot(outlier.shape=23,
               outlier.size=1) +
  facet_wrap(~p, ncol=3,labeller="label_both")+
  xlab("Sample Size")+
  ylab(expression('Correctly recovered ancestors \nover total causal ancestors'))+
  theme(plot.margin = unit(c(0,0,0,5), "mm"))


Plot <- ggarrange(ANC1, ANC2, nrow=2, common.legend = TRUE, legend="right")

ggsave(
  "Plots/CAMvsTree_Singlerooted_Ancestors.pdf",
  plot = Plot,
  device = "pdf",
  path = NULL,
  scale = 1,
  width = 210,
  height = 140,
  units = c("mm"),
  dpi = 500,
  limitsize = TRUE
)

