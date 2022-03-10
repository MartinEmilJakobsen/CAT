#Tree generation functions
Generate_tree_adjacencymatrix <- function(type,p){
  Adjacency <- matrix(0,ncol=p,nrow=p)
  if(type == 1){
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
  }
  else if(type ==2){
        for(j in 2:p){
          i <- sample(seq(1,j-1,1),size=1)
          Adjacency[i,j] <- 1
          
        }
    
  } else{
    stop("type not recognized")
  }
  
  return(Adjacency)
}

#Take Adjacency matrix and create adjacency data-format with edge rows parent,node 
Generate_tree_adjacency_dat <- function(Adjacency,p){
  
  roots <- which(colSums(Adjacency)==0)
  
  Adj.Dat <- data.frame(parent=0,node = roots) 
  
  for(i in 1:p){
    for(j in 1:p){
      if(Adjacency[i,j]!=0)
      {
        Adj.Dat <- rbind(Adj.Dat,data.frame(parent=i,node = j)  )
      }
    }}
  Adj.Dat <- Adj.Dat %>% arrange(node)
  return(Adj.Dat)
}

#Calculate ancestor lists from Adjacency Matrix
Find_Ancestors <- function(Adjacency,p){ 
  expAdjMat <- {expm(Adjacency)} %>% as.matrix()
  #expAdjMat[lower.tri(expAdjMat,diag=TRUE)] <- 0
  expAdjMat[diag(p)==1] <- 0
  
  An <- list()
  for(i in 1:p){
    An[[i]] <-  which(expAdjMat[,i]!=0)
  }
  return(tibble(node=1:p,Ancestors=An))
}


#Generate noise innovations for simulation

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

#Generate GP causal function SCM data

Generate_GP_Data <- function(noisedist,noisevariances,samplesize,Adj.Dat,bandwidth){
  Noise <- Generate_noise(samplesize,noisedist,noisevariances,p)
  l <- bandwidth 
  Data <- matrix(0,ncol=p,nrow=samplesize) 
  DatawoNoise <- matrix(0,ncol=p,nrow=samplesize) 
  Data[,1] <- Noise[,1]
  for(nod in 2:p){
    d = abs(outer(Data[,Adj.Dat[nod,1]],Data[,Adj.Dat[nod,1]],"-")) # compute distance matrix, d_{ij} = |x_i - x_j|
    Sigma_SE = exp(-d^2/(2*l^2)) # squared exponential kernel
    rm(d)
    Data[,nod] <- mgcv:::rmvn(1, rep(0,samplesize), Sigma_SE) + Noise[,nod]
    rm(Sigma_SE)
  }
  colnames(Data) <- paste0("V",seq(1,p,1),sep="")
  Data <- Data %>% as_tibble()
  return(Data)
}



#Gen CausalFunctions as GAM models.

Generate_CausalFunctions <- function(Adj.Dat,bandwidth){
  
  CausalFunctions <- Adj.Dat %>%  as_tibble() %>% 
    mutate(CausalFunction = list(list()))
  
  for( i in 1:nrow(CausalFunctions)){
    from <- CausalFunctions[i,"parent"]
    to <- CausalFunctions[i,"node"]
    from_name <- paste0("V",from)
    to_name <- paste0("V",to)
    TempData <- tibble(!!from_name:=seq(-5,5,10/100),!!to_name:=0)
    d = abs(outer(TempData[,1] %>% pull,TempData[,1] %>% pull,"-")) # compute distance matrix, d_{ij} = |x_i - x_j|
    Sigma_SE = exp(-d^2/(2*1^2)) # squared exponential kernel
    rm(d)
    TempData[,2] <- mgcv:::rmvn(1, rep(0,nrow(TempData)), Sigma_SE)
    form <- formula(paste0(to_name,"~s(",from_name,",bs='tp')"))
    CausalFunctions[i,"CausalFunction"][[1]] <- list(gam(form,data=TempData))
    
  }
  
  return(CausalFunctions)
}
#Gen Data from Gam Causal Functions

Generate_Data_From_Causal_Functions <- function(noisedist,noisevariances,samplesize,Adj.Dat,bandwidth,p,CausalFunctions){
  
  l <- bandwidth
  Data <- Generate_noise(samplesize,noisedist,noisevariances,p) 
  colnames(Data) <- paste0("V",seq(1,p,1),sep="")
  Data <- Data %>% as_tibble()
  
  for(nod in 2:p){
    parents <-Adj.Dat %>% filter(node==nod) %>% select(parent) %>%  pull
    for (pare in parents){
      Data[,nod] <- Data[,nod] + predict(CausalFunctions %>% filter(parent==pare,node==nod) %>% select(CausalFunction) %>% pull %>% {.[[1]]},newdata=Data )
    }
  }
  return(Data)
}




CAT <- function(dat,noise,workers,numbasisfnct = NULL)
 {
  
  #Saving original column names and setting standard column names for dat-processing
  colNames <- names(dat)
  names(dat) <- paste0("X",seq(1,ncol(dat),1)) 
  p <- ncol(dat)
  
  if(is.null(numbasisfnct)){
               if(noise =="G"){
               f.edgeweight <- function(from,to){
                 form <- formula(paste0("X",to,"~","s(X",from,",bs='tp')"))
                 Cond_exp <- gam(form, data = dat)
                 form <- paste0("dat$X",to,"- predict(Cond_exp,newdata=dat)")
                 Residual <- eval(parse(text=form))
                 var(Residual)
                 } 
               } else if(noise=="NA"){  
               f.edgeweight <- function(from,to){
                 form <- formula(paste0("X",to,"~","s(X",from,",bs='tp')"))
                 Cond_exp <- gam(form, data = dat)
                 form <- paste0("dat$X",to,"- predict(Cond_exp,newdata=dat)")
                 Residual <- eval(parse(text=form))
                 KLentropy(unique(Residual),k=4)$Estimate
                 } 
               } else { 
                 stop(paste0("Error: Expected noise='G'/'NA' but recieved noise='",noise,"'"))
               }
  } else {
                if(noise =="G"){
                  f.edgeweight <- function(from,to){
                    form <- formula(paste0("X",to,"~","s(X",from,", k=numbasisfnct, bs='tp')"))
                    Cond_exp <- gam(form, data = dat)
                    form <- paste0("dat$X",to,"- predict(Cond_exp,newdata=dat)")
                    Residual <- eval(parse(text=form))
                    var(Residual)
                  } 
                } else if(noise=="NA"){  
                  f.edgeweight <- function(from,to){
                    form <- formula(paste0("X",to,"~","s(X",from,", k=numbasisfnct, bs='tp')"))
                    Cond_exp <- gam(form, data = dat)
                    form <- paste0("dat$X",to,"- predict(Cond_exp,newdata=dat)")
                    Residual <- eval(parse(text=form))
                    KLentropy(unique(Residual),k=4)$Estimate
                  } 
                } else { 
                  stop(paste0("Error: Expected noise='G'/'NA' but recieved noise='",noise,"'"))
                }    
  }
    
    
  if(noise =="G"){
    f.targetweight <- function(node){
      varstring <- paste0("X",node)
      variable <- dat %>% select(all_of(varstring)) %>% pull
      var(variable)
      
    } 
  } else if(noise=="NA"){  
    f.targetweight <- function(node){
      varstring <- paste0("X",node)
      variable <- dat %>% select(all_of(varstring)) %>% pull
      KLentropy(variable,k=4)$Estimate
    } 
  }
  
  
  #All node combinations
  Edges <- expand_grid(from=seq(1,p,1),to=seq(1,p,1)) %>% filter(from != to)
  Nodes <- expand_grid(node=seq(1,p,1))
  
  #Calculating the Edge Scores
  options(future.rng.onMisuse="ignore") #Ignore seed messages
  
  #message("\n Computing EdgeWeights:")
  if(workers>1){
  plan(multisession, workers = workers)
  EdgeWeights <- Edges %>% 
    mutate(EdgeWeight = future_pmap_dbl(.l=list(from,to),
                                 .f=f.edgeweight,.progress=TRUE)) 
  future:::ClusterRegistry("stop")
  }  else{
  EdgeWeights <- Edges %>% 
    mutate(EdgeWeight = pmap_dbl(.l=list(from,to),
                                          .f=f.edgeweight)) 
  }
  
  #message("\n Computing TargetWeights:")
  if(workers>1){
  plan(multisession, workers = workers)
  TargetWeights <- Nodes %>% 
    mutate(TargetWeight = future_pmap_dbl(.l=list(node),
                                        .f=f.targetweight,.progress=TRUE)) 
  future:::ClusterRegistry("stop")
  } else{
    TargetWeights <- Nodes %>% 
      mutate(TargetWeight = pmap_dbl(.l=list(node),
                                            .f=f.targetweight)) 
  }
  
  
  if(noise=="G"){
    EdgeWeightMatrix <- left_join(EdgeWeights,TargetWeights,by=c("to"="node")) %>% 
      mutate(NegEdgeScore = -log(EdgeWeight/TargetWeight) ,
             NegEdgeScore = NegEdgeScore+2*abs(min(NegEdgeScore))
             ) %>% 
      select(from,to,NegEdgeScore) %>% 
      spread(key=to,value=NegEdgeScore,fill = 0) %>% 
      arrange(from) %>% #increasing order, so first row is node 1
      select(-from) %>% 
      as.matrix
  } else if(noise =="NA"){
    EdgeWeightMatrix <- left_join(EdgeWeights,TargetWeights,by=c("to"="node")) %>% 
      mutate(NegEdgeScore = TargetWeight-EdgeWeight,
             NegEdgeScore = NegEdgeScore+2*abs(min(NegEdgeScore))) %>% 
      select(from,to,NegEdgeScore) %>% 
      spread(key=to,value=NegEdgeScore,fill = 0) %>% 
      arrange(from) %>% #increasing order, so first row is node 1
      select(-from) %>% 
      as.matrix
    
  }
 
  EdgeScoresAsGraphNEL<-as(EdgeWeightMatrix, "graphNEL")
  Optimal<-edmondsOptimumBranching(EdgeScoresAsGraphNEL)
  

  
  AdjMatrix <- matrix(rep(0,ncol(dat)*ncol(dat)),ncol=ncol(dat))
  
  for(i in 1:ncol(Optimal$edgeList)){
    from <- Optimal$edgeList["from",i] %>% as.integer
    to <- Optimal$edgeList["to",i] %>% as.integer
    AdjMatrix[from,to] <- 1
  }
  rownames(EdgeWeightMatrix)<-colNames
  colnames(EdgeWeightMatrix)<-colNames
  
  AdjMatrix <- as(AdjMatrix,"dgCMatrix")
  
  reslist <- list(EdgeWeightMatrix=EdgeWeightMatrix,AdjMatrix=AdjMatrix,EdgeWeights=EdgeWeights,TargetWeights=TargetWeights,InputDatagraphNEL = EdgeScoresAsGraphNEL, Optimal =Optimal )
  
  return(reslist)
}




#### Hypothesis_Test Takes multiple simple hypothesis and test each single one, and returns a vector of TRUE/FALSE 
HypothesisTest <- function(dat,hyp,alpha,numbasisfnct = NULL){
  colNames <- names(dat)
  names(dat) <- paste0("X",seq(1,ncol(dat),1)) 
  p <- ncol(dat)
  Data1 <- dat[1:floor(nrow(dat)/2),]
  Data2 <- dat[(floor((nrow(dat)/2)+1)):nrow(dat),]
  n <- nrow(Data1)
  quantile = qnorm(alpha/(2*(p*(p-1))),
                   mean=0,
                   sd= 1,lower.tail = FALSE)
  CI <- expand_grid(from=seq(1,p,1),to=seq(1,p,1)) %>%  
    filter(from!= to) %>% 
    rowwise() %>% 
    mutate(M = list(pull((Data1[,to]-predict(gam(formula(paste0("X",to,"~","s(X",from,",bs='tp',k=ifelse(is.null(numbasisfnct),-1,numbasisfnct))")), data = Data2),
                                             newdata=Data1[,paste0("X",from)]))^2)),
           V = list(pull((Data1[,to]-mean(pull(Data1[,to])))^2))) %>% 
    unnest(cols=c("M","V")) %>% 
    ungroup() %>% 
    group_by(from,to) %>% 
    summarize(mu = mean(M),
              Sigma_M = mean(M^2-mu^2),
              nu = mean(V),
              Sigma_V = mean(V^2-nu^2),
              Sigma_MV = mean(M*V-mu*nu),
              sigma = Sigma_M/mu^2 + Sigma_V/nu^2 - 2*Sigma_MV/(mu*nu),
              lower = (1/2)*log(mu/nu) - quantile*sigma/(2*sqrt(n)),
              upper = (1/2)*log(mu/nu) + quantile*sigma/(2*sqrt(n)),
              edgeweight = (1/2)*log(mu/nu),
              add = quantile*sigma/(2*sqrt(n)),
              .groups = 'drop') %>% 
    mutate(lowerPOS = lower + 2*abs(min(lower)),
           upperPOS = upper + 2*abs(min(lower)))


  LowerWeigthedAdjacency_Bonferroni <- matrix(rep(NA,p*p),ncol=p)
  UpperWeigthedAdjacency_Bonferroni <- matrix(rep(NA,p*p),ncol=p)
   
  for(i in 1:p){
    for (j in 1:p){
      if (i != j){
        LowerWeigthedAdjacency_Bonferroni[j,i] <- CI %>% 
          filter(from == j,to == i) %>% 
          select(lowerPOS) %>% 
          pull()
        UpperWeigthedAdjacency_Bonferroni[j,i] <- CI %>% 
          filter(from == j,to == i) %>% 
          select(upperPOS) %>% 
          pull()
      }
      else {
        LowerWeigthedAdjacency_Bonferroni[j,i] <- 0
        UpperWeigthedAdjacency_Bonferroni[j,i] <- 0
      }
    }
  }

  LowerEdInput_Bonferroni <- paste(c(t(LowerWeigthedAdjacency_Bonferroni)), collapse = " ")
  UpperEdInput_Bonferroni <- paste(c(t(UpperWeigthedAdjacency_Bonferroni)), collapse = " ")

  HypInput <- paste(c(t(hyp)), collapse = " ")

  PythonInput_Bonferroni = paste(p,LowerEdInput_Bonferroni,UpperEdInput_Bonferroni,HypInput)

  Test_Exact_Bonferroni = system(paste('python3 PythonMultipleTests_ExactQuery.py',PythonInput_Bonferroni), intern = TRUE) %>%
    strsplit(., split = " ") %>%
    unlist() %>%
    as.integer() %>%
    as.vector()

  Test_Asym_Bonferroni = system(paste('python3 PythonMultipleTests_AsymQuery.py',PythonInput_Bonferroni), intern = TRUE) %>% 
    strsplit(., split = " ") %>% 
    unlist() %>% 
    as.integer() %>% 
    as.vector()
  

  hyp <- hyp %>% 
     mutate( Test_Exact_Bonferroni = Test_Exact_Bonferroni,
             Test_Asym_Bonferroni = Test_Asym_Bonferroni
           )
  
  return(hyp)
}




