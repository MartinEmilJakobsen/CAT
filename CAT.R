
CAT <- function(dat,noise="G",workers=1,numbasisfnct = NULL, pvalCutoff = 1, crossfit = FALSE)
{
  colNames <- names(dat)
  names(dat) <- paste0("X",seq(1,ncol(dat),1)) 
  p <- ncol(dat)
  
  CrossFit <- function(DataTrainCondExp,DataTrainEdgeWeight,noise,workers,numbasisfnct){
  #Saving original column names and setting standard column names for dat-processing
  if(is.null(numbasisfnct)){
    if(noise =="G"){
      f.edgeweight <- function(from,to){
        form <- formula(paste0("X",to,"~","s(X",from,",bs='tp')"))
        Cond_exp <- gam(form, data = DataTrainCondExp)
        form <- paste0("DataTrainEdgeWeight$X",to,"- predict(Cond_exp,newdata=DataTrainEdgeWeight)")
        Residual <- eval(parse(text=form))
        data.frame(pval = summary.gam(Cond_exp)$s.pv,EdgeWeight = var(Residual))
      } 
    } else if(noise=="NA"){  
      f.edgeweight <- function(from,to){
        form <- formula(paste0("X",to,"~","s(X",from,",bs='tp')"))
        Cond_exp <- gam(form, data = DataTrainCondExp)
        form <- paste0("DataTrainEdgeWeight$X",to,"- predict(Cond_exp,newdata=DataTrainEdgeWeight)")
        Residual <- eval(parse(text=form))
        data.frame(pval = summary.gam(Cond_exp)$s.pv,EdgeWeight = KLentropy(unique(Residual),k=4)$Estimate)
      } 
    } else { 
      stop(paste0("Error: Expected noise='G'/'NA' but recieved noise='",noise,"'"))
    }
  } else {
    if(noise =="G"){
      f.edgeweight <- function(from,to){
        form <- formula(paste0("X",to,"~","s(X",from,", k=numbasisfnct, bs='tp')"))
        Cond_exp <- gam(form, data = DataTrainCondExp)
        form <- paste0("DataTrainEdgeWeight$X",to,"- predict(Cond_exp,newdata=DataTrainEdgeWeight)")
        Residual <- eval(parse(text=form))
        data.frame(pval = summary.gam(Cond_exp)$s.pv,EdgeWeight = var(Residual))
      } 
    } else if(noise=="NA"){  
      f.edgeweight <- function(from,to){
        form <- formula(paste0("X",to,"~","s(X",from,", k=numbasisfnct, bs='tp')"))
        Cond_exp <- gam(form, data = DataTrainCondExp)
        form <- paste0("DataTrainEdgeWeight$X",to,"- predict(Cond_exp,newdata=DataTrainEdgeWeight)")
        Residual <- eval(parse(text=form))
        data.frame(pval = summary.gam(Cond_exp)$s.pv,EdgeWeight = KLentropy(unique(Residual),k=4)$Estimate)
      } 
    } else { 
      stop(paste0("Error: Expected noise='G'/'NA' but recieved noise='",noise,"'"))
    }    
  }
  
  
  if(noise =="G"){
    f.targetweight <- function(node){
      varstring <- paste0("X",node)
      variable <- DataTrainEdgeWeight %>% select(all_of(varstring)) %>% pull
      var(variable)
      
    } 
  } else if(noise=="NA"){  
    f.targetweight <- function(node){
      varstring <- paste0("X",node)
      variable <- DataTrainEdgeWeight %>% select(all_of(varstring)) %>% pull
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
      mutate(res = future_pmap(.l=list(from,to),
                                          .f=f.edgeweight,.progress=TRUE)) 
    future:::ClusterRegistry("stop")
  }  else{
    EdgeWeights <- Edges %>% 
      mutate(res = pmap(.l=list(from,to),
                                   .f=f.edgeweight)) 
  }
  EdgeWeights <- EdgeWeights %>%  unnest(cols=c("res"))
  
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
  
  return(left_join(EdgeWeights,TargetWeights,by=c("to"="node")))
  
  
  
  }
  
  if(crossfit == TRUE){
    idx <- sample(nrow(dat),nrow(dat)/2,replace=FALSE)
    dat1 <- dat[idx,]
    dat2 <- dat[-idx,]
    
    ET.Weights1 <- CrossFit(DataTrainCondExp=dat1,DataTrainEdgeWeight=dat2,noise=noise,workers=workers,numbasisfnct=numbasisfnct)
    ET.Weights2 <- CrossFit(DataTrainCondExp=dat2,DataTrainEdgeWeight=dat1,noise=noise,workers=workers,numbasisfnct=numbasisfnct)
    
    ET.Weights <- bind_rows(ET.Weights1,ET.Weights2) %>% 
      group_by(from,to) %>% 
      summarize_all(.funs=mean)
    
  } else {
    ET.Weights <- CrossFit(DataTrainCondExp=dat,DataTrainEdgeWeight=dat,noise=noise,workers=workers,numbasisfnct=numbasisfnct)
  }
  
  
  if(noise=="G"){
    EdgeWeightMatrix <- ET.Weights %>% 
      ungroup() %>% 
      mutate(NegEdgeScore = -log(EdgeWeight/TargetWeight) ,
             NegEdgeScore = NegEdgeScore+2*abs(min(NegEdgeScore))
      ) %>% 
      select(from,to,NegEdgeScore) %>% 
      spread(key=to,value=NegEdgeScore,fill = 0) %>% 
      arrange(from) %>% #increasing order, so first row is node 1
      select(-from) %>% 
      as.matrix
  } else if(noise =="NA"){
    EdgeWeightMatrix <- ET.Weights %>% 
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
  
  prunelist <- ET.Weights %>% filter(pval>pvalCutoff)
  
  for (i in 1:nrow(prunelist)){
    AdjMatrix[prunelist[i,"from"] %>% as.integer() ,prunelist[i,"to"] %>%  as.integer()] <- 0
  }
  
  AdjMatrix <- as(AdjMatrix,"dgCMatrix")

  reslist <- list(EdgeWeightMatrix=EdgeWeightMatrix,AdjMatrix=AdjMatrix, Optimal =Optimal )
  
  return(reslist)
}

HypothesisTest <- function(dat,hypInput,level,queryscheme,numbasisfnct = NULL){
  alpha <- level
  colNames <- names(dat)
  names(dat) <- paste0("X",seq(1,ncol(dat),1)) 
  
  #translating hypothesis to new names
  translation <- data.frame(original = as.character(colNames),new=seq(1,ncol(dat),1))
  hypInput <- hypInput %>% mutate(from=as.character(from),to=as.character(to)) %>% 
    rowwise() %>% mutate(from = translation[which(translation[,1]==from),2],
                         to = translation[which(translation[,1]==to),2])

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
  
  HypInputPython <- paste(c(t(hypInput)), collapse = " ")
  
  PythonInput_Bonferroni = paste(queryscheme,p,LowerEdInput_Bonferroni,UpperEdInput_Bonferroni,HypInputPython)  

  
  Test = system(paste('python3 PythonTest.py',PythonInput_Bonferroni), intern = TRUE) %>%
    strsplit(., split = " ") %>%
    unlist() %>%
    as.integer() %>%
    as.vector()
  
  return(Test)
}
