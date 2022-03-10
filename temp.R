p <- 16
qmvnorm(0.95,sigma=diag(rep(1,p*p-1)),tail=c("lower.tail"))
qmvnorm(0.95,sigma=diag(rep(1,p*p-1)),tail=c("upper.tail"))
qmvnorm(0.95,sigma=diag(rep(1,p*p-1)),tail=c("both.tails"))

qnorm((1-0.05/(2*p*(p-1))),mean=0,sd=1)
      
  
quantile = qnorm(alpha/(2*(p*(p-1))),
                 mean=0,
                 sd= 1,lower.tail = FALSE)


CI <- expand_grid(from=seq(1,p,1),to=seq(1,p,1)) %>%  
  filter(from!= to) %>% 
  rowwise() %>% 
  mutate(M = list(pull((Data1[,to]-predict(gam(formula(paste0("X",to,"~","s(X",from,",bs='tp')")), data = Data2),
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


#New CI


M <- expand_grid(to=seq(1,p,1),from=seq(1,p,1)) %>%  
  filter(from!= to) %>% 
  rowwise() %>% 
  mutate(M = list(pull((Data1[,to]-predict(gam(formula(paste0("X",to,"~","s(X",from,",bs='tp')")), data = Data2),
                                           newdata=Data1[,paste0("X",from)]))^2))) %>% 
  unnest(cols=c("M"))

mu <- M %>% group_by(from,to) %>% summarize(mu = mean(M)) %>% ungroup()


V <- expand_grid(to=seq(1,p,1)) %>%  
  rowwise() %>% 
  mutate(V = list(pull((Data1[,to]-mean(pull(Data1[,to])))^2))) %>% 
  unnest(cols=c("V"))

nu <- V %>% group_by(to) %>% summarize(nu = mean(V)) %>% ungroup()


# in matrix norm with base R covariance matrix

M_mat <- M %>%  
  mutate(var_name = paste0("M",from,to)) %>%
  mutate(rep = rep(seq(1,n,1),p*(p-1))) %>% 
  select(var_name,M,rep) %>% 
  pivot_wider(names_from= var_name, values_from =M)


V_mat <- V %>%  
  mutate(var_name = paste0("V",to)) %>%
  mutate(rep = rep(seq(1,n,1),p)) %>% 
  select(var_name,V,rep) %>% 
  pivot_wider(names_from= var_name, values_from =V)

designmatrix <- left_join(M_mat,V_mat,by=c("rep")) %>% select(-rep)

Sigma <- cov(designmatrix)

Dg <- matrix(rep(NA,p*(p-1)*(p*(p-1)+p)),ncol=p*(p-1)+p)

r <- 1
column <- 1
for (to1 in seq(1,p,1)){
  for (from1 in seq(1,p,1)){
    column <- 1
    if (from1 != to1){
      #deriv row
      for (to2 in seq(1,p,1)){
        for (from2 in seq(1,p,1)){
          if( from2 != to2){
            if (from1==from2 & to1 == to2){
              Dg[r,column] <- 1/ (2*(mu %>%  filter(from == from1, to == to1) %>% select(mu) %>%  pull))
              column <- column + 1
            }else {
              Dg[r,column] <- 0
              column <- column + 1
          } 
          }
        }
      }
      for (to3 in seq(1,p,1)){
        if(to3 == to1){
          Dg[r,column] <- -1/(2*(nu %>%  filter(to == to1) %>% select(nu) %>%  pull))
        } else{
          Dg[r,column] <- 0
        }
        column = column + 1
      }
      r <- r +1
    }
  }
}

q <- qmvnorm(0.95,sigma=diag(rep(1,p*p-1)),tail=c("both.tails"))

q_alpha <- rep(q$quantile,p*(p-1))

Simult.Data <- expand_grid(to=seq(1,p,1),from=seq(1,p,1)) %>% 
  filter(from != to) %>% 
  select(from,to) %>% 
  mutate(add.simult = sqrtm(Dg %*% Sigma %*% t(Dg))[,1] %*% q_alpha / sqrt(n) %>%  unlist %>% as.vector) 



CI <- left_join(CI,Simult.Data,cols=c("from","to")) %>% 
  mutate(lowerPOS.simlt = edgeweight - add.simult + 2*abs(min(edgeweight - add.simult)),
         upperPOS.simlt = edgeweight + add.simult + 2*abs(min(edgeweight - add.simult))) %>% 
  mutate(lower_diff = lowerPOS.simlt-lowerPOS) 

#### FROM ANALYSISOFDATA TESTING

SummaryDat <- readRDS(file="Data/Data_Testing_Multivariate_2201061730.RDS")

SummaryDat %>% unnest(cols=c("Test")) %>% 
  mutate(TestVal = ifelse(Test==TRUE, 1, 0)) %>% 
  group_by(N,p,TreeType,HypothesisType,HypothesisLocation) %>% 
  summarise( from=mean(from),to=mean(to),val = mean(TestVal)) %>% 
  mutate(Power = ifelse(HypothesisType=="Missing",1-val,NA),
         Level = ifelse(HypothesisType=="Present",val,NA)) %>% 
  print(n=100)


SummaryDat %>% unnest(cols=c("Test")) %>% filter(HypothesisType == "Missing")%>% print(n=500)






# 8 Nodes 50k samples, Location="Root" only
SummaryDat <- readRDS(file="Data/Data_Testing_Multivariate_50k_2201071428.RDS")

SummaryDat %>% unnest(cols=c("Test")) %>% 
  mutate(Reject = ifelse(Test==TRUE, 0, 1)) %>% 
  group_by(N,p,TreeType,HypothesisType,HypothesisLocation) %>% 
  summarise( from=mean(from),to=mean(to),Stat_Val = mean(Reject)) %>%  
  mutate(Stat = ifelse(HypothesisType == "Missing","Power","Level")) %>% 
  arrange(HypothesisType,HypothesisLocation) %>%  
  select(N,p,TreeType,HypothesisType,HypothesisLocation,from,to,Stat,Stat_Val) %>% 
  print(n=100)

# 8 nodes


SummaryDat <- readRDS(file="Data/Data_Testing_Multivariate_2201071650.RDS")

SummaryDat %>% unnest(cols=c("Test")) %>% 
  mutate(Reject = ifelse(Test==TRUE, 0, 1)) %>% 
  group_by(N,p,TreeType,HypothesisType,HypothesisLocation) %>% 
  summarise( from=mean(from),to=mean(to),Stat_Val = mean(Reject)) %>%  
  mutate(Stat = ifelse(HypothesisType == "Missing","Power","Level")) %>% 
  arrange(HypothesisType,HypothesisLocation) %>%  
  select(N,p,TreeType,HypothesisType,HypothesisLocation,from,to,Stat,Stat_Val) %>% 
  print(n=100)


# 8 nodes after fix (still wrong - forgot to update functions.R)



SummaryDat <- readRDS(file="Data/Data_Testing_Multivariate_2201081455.RDS")

SummaryDat %>% unnest(cols=c("Test")) %>% 
  mutate(Reject = ifelse(Test==TRUE, 0, 1)) %>% 
  group_by(N,p,TreeType,HypothesisType,HypothesisLocation) %>% 
  summarise( from=mean(from),to=mean(to),Stat_Val = mean(Reject)) %>%  
  mutate(Stat = ifelse(HypothesisType == "Missing","Power","Level")) %>% 
  arrange(HypothesisType,HypothesisLocation) %>%  
  select(N,p,TreeType,HypothesisType,HypothesisLocation,from,to,Stat,Stat_Val) %>% 
  print(n=100)


# 8 nodes real fix



SummaryDat <- readRDS(file="Data/Data_Testing_Multivariate_2201081702.RDS")

SummaryDat %>% unnest(cols=c("Test")) %>% 
  mutate(Reject = ifelse(Test==TRUE, 0, 1)) %>% 
  group_by(N,p,TreeType,HypothesisType,HypothesisLocation) %>% 
  summarise( from=mean(from),to=mean(to),Stat_Val = mean(Reject)) %>%  
  mutate(Stat = ifelse(HypothesisType == "Missing","Power","Level")) %>% 
  arrange(HypothesisType,HypothesisLocation) %>%  
  select(N,p,TreeType,HypothesisType,HypothesisLocation,from,to,Stat,Stat_Val) %>% 
  print(n=100)




# 8 nodes v2

SummaryDat <- readRDS(file="Data/Data_Testing_Multivariate_v2_2201091544.RDS")

SummaryDat %>% unnest(cols=c("Test")) %>% 
  mutate(Reject = ifelse(Test==TRUE, 0, 1)) %>% 
  group_by(N,p,TreeType,HypothesisType,HypothesisLocation) %>% 
  summarise( from=mean(from),to=mean(to),Stat_Val = mean(Reject)) %>%  
  mutate(Stat = ifelse(HypothesisType == "Missing","Power","Level")) %>% 
  arrange(HypothesisType,HypothesisLocation) %>%  
  select(N,p,TreeType,HypothesisType,HypothesisLocation,from,to,Stat,Stat_Val) %>% 
  print(n=100)


# 8 nodes v2 new

SummaryDat <- readRDS(file="Data/Data_Testing_Multivariate_v2_2201101622.RDS")

SummaryDat %>% unnest(cols=c("Test")) %>% 
  mutate(Reject = ifelse(Test==TRUE, 0, 1)) %>% 
  mutate(Stat = ifelse(HypothesisType == "Missing","Power","Level")) %>% 
  group_by(N,p,TreeType,HypothesisType,HypothesisLocation,Stat) %>% 
  summarise( from=mean(from),to=mean(to),Stat_Val = mean(Reject)) %>%  
  arrange(HypothesisType,HypothesisLocation) %>%  
  select(N,p,TreeType,HypothesisType,HypothesisLocation,from,to,Stat,Stat_Val) %>% 
  print(n=100)

SummaryDat %>%  filter(N == 25000) %>% 
  unnest(cols=c("Test")) %>% 
  mutate(Stat = ifelse( (IsEdgePresent == "TRUE" & Hyp == 1) | (IsEdgePresent == FALSE & Hyp == -1 ),"Level","Power")) %>% 
  group_by(N,p,TreeType,IsEdgePresent,dist,Hyp,Stat) %>% 
  summarise( from=mean(from),to=mean(to),Stat_Val = mean(Test), n = n()) %>%  
  arrange(N,Stat,dist) %>% 
  print(n=200)


##### 8 NODES LARGE SAMPLESIZE

SummaryDat <- readRDS(file="Data/Data_Testing_Multivariate_v2_2201101733.RDS")


dat1 <- SummaryDat %>%  # filter(N == 25000) %>% 
  unnest(cols=c("Test")) %>% 
  mutate(Stat = ifelse( (IsEdgePresent == "TRUE" & Hyp == 1) | (IsEdgePresent == FALSE & Hyp == -1 ),"Level","Power")) %>% 
  mutate(dist = case_when(dist < 0 ~ "<0", dist == 0 ~ "0", dist > 0 ~ ">0")) %>% 
  group_by(N,p,TreeType,IsEdgePresent,dist,Hyp,Stat) %>% 
  summarise( from=mean(from),
             to=mean(to),
             Stat_Val = mean(Test) , 
             n = n(), 
             lower = Stat_Val-qnorm(0.975,mean=0,sd=1)*sqrt(Stat_Val*(1-Stat_Val)/n),
             upper = Stat_Val+qnorm(0.975,mean=0,sd=1)*sqrt(Stat_Val*(1-Stat_Val)/n)) %>%  
  arrange(IsEdgePresent,N,Stat,dist) %>% 
  print(n=200)

dat1 %>%  filter(IsEdgePresent == FALSE, Stat == "Power") %>% 
  mutate(dist= factor(dist)) %>% 
  ggplot(data=.)+
  geom_line(aes(x=N,y=Stat_Val,group = dist, color = dist))+
  ggtitle(expression("Power -- Truth: NOT j" %->% "i,   H_0: j" %->% "i"))+
  geom_ribbon(aes(x=N,ymin=lower,ymax=upper,group = dist),alpha=0.3)


dat1 %>%  filter(IsEdgePresent == TRUE, Stat == "Power") %>% 
  mutate(dist= factor(dist)) %>% 
  ggplot(data=.)+
  geom_line(aes(x=N,y=Stat_Val,group = dist, color = dist))+
  ggtitle(expression("Power -- Truth: j" %->% "i,   H_0: NOT j" %->% "i"))



# 8 Nodes even more samplesize
SummaryDat <- readRDS(file="Data/Data_Testing_Multivariate_v2_2201111146.RDS")


dat <- SummaryDat %>%  # filter(N == 25000) %>% 
  unnest(cols=c("Test")) %>% 
  mutate(Stat = ifelse( (IsEdgePresent == "TRUE" & Hyp == 1) | (IsEdgePresent == FALSE & Hyp == -1 ),"Level","Power")) %>% 
  mutate(dist = case_when(dist < 0 ~ "<0", dist == 0 ~ "0", dist > 0 ~ ">0")) %>% 
  group_by(N,p,TreeType,IsEdgePresent,dist,Hyp,Stat) %>% 
  summarise( from=mean(from),to=mean(to),Stat_Val = mean(Test), n = n()) %>%  
  arrange(IsEdgePresent,N,Stat,dist) %>% 
  print(n=200)

dat %>%  filter(IsEdgePresent == FALSE, Stat == "Power") %>% 
  mutate(dist= factor(dist)) %>% 
  ggplot(data=.)+
  geom_line(aes(x=N,y=Stat_Val,group = dist, color = dist))+
  ggtitle(expression("Power -- Truth: NOT j" %->% "i,   H_0: j" %->% "i"))


dat %>%  filter(IsEdgePresent == TRUE, Stat == "Power") %>% 
  mutate(dist= factor(dist)) %>% 
  ggplot(data=.)+
  geom_line(aes(x=N,y=Stat_Val,group = dist, color = dist))+
  ggtitle(expression("Power -- Truth: j" %->% "i,   H_0: NOT j" %->% "i"))



#8 nodes half variance


SummaryDat <- readRDS(file="Data/Data_Testing_Multivariate_v2_2201112257.RDS")


dat <- SummaryDat %>%  # filter(N == 25000) %>% 
  unnest(cols=c("Test")) %>% 
  mutate(Stat = ifelse( (IsEdgePresent == "TRUE" & Hyp == 1) | (IsEdgePresent == FALSE & Hyp == -1 ),"Level","Power")) %>% 
  mutate(dist = case_when(dist < 0 ~ "<0", dist == 0 ~ "0", dist > 0 ~ ">0")) %>% 
  group_by(N,p,TreeType,IsEdgePresent,dist,Hyp,Stat) %>% 
  summarise( from=mean(from),to=mean(to),Stat_Val = mean(Test), n = n()) %>%  
  arrange(IsEdgePresent,N,Stat,dist) %>% 
  print(n=200)

dat %>%  filter(IsEdgePresent == FALSE, Stat == "Power") %>% 
  mutate(dist= factor(dist)) %>% 
  ggplot(data=.)+
  geom_line(aes(x=N,y=Stat_Val,group = dist, color = dist))+
  ggtitle(expression("Power -- Truth: NOT j" %->% "i,   H_0: j" %->% "i"))


dat %>%  filter(IsEdgePresent == TRUE, Stat == "Power") %>% 
  mutate(dist= factor(dist)) %>% 
  ggplot(data=.)+
  geom_line(aes(x=N,y=Stat_Val,group = dist, color = dist))+
  ggtitle(expression("Power -- Truth: j" %->% "i,   H_0: NOT j" %->% "i"))


#full and halfvariance
rbind(dat1 %>% mutate(type= "Half var"),dat %>% mutate(type="Norm var")) %>% 
  filter(N <= 100000) %>% 
  filter(IsEdgePresent == FALSE, Stat == "Power") %>% 
  mutate(dist= factor(dist)) %>% 
  ggplot(data=.)+
  geom_line(aes(x=N,y=Stat_Val,group = interaction(dist,type), color = dist,linetype =type))+
  ggtitle(expression("Power -- Truth: NOT j" %->% "i,   H_0: j" %->% "i"))



# 16 nodes

SummaryDat <- readRDS(file="Data/Data_Testing_Multivariate_v2_2201121134.RDS")


dat <- SummaryDat %>%  # filter(N == 25000) %>% 
  unnest(cols=c("Test")) %>% 
  mutate(Stat = ifelse( (IsEdgePresent == "TRUE" & Hyp == 1) | (IsEdgePresent == FALSE & Hyp == -1 ),"Level","Power")) %>% 
  mutate(dist = case_when(dist < 0 ~ "<0", dist == 0 ~ "0", dist > 0 ~ ">0")) %>% 
  group_by(N,p,TreeType,IsEdgePresent,dist,Hyp,Stat) %>% 
  summarise( from=mean(from),to=mean(to),Stat_Val = mean(Test), n = n()) %>%  
  arrange(IsEdgePresent,N,Stat,dist) %>% 
  print(n=200)

dat %>%  filter(IsEdgePresent == FALSE, Stat == "Power") %>% 
  mutate(dist= factor(dist)) %>% 
  ggplot(data=.)+
  geom_line(aes(x=N,y=Stat_Val,group = dist, color = dist))+
  ggtitle(expression("Power -- Truth: NOT j" %->% "i,   H_0: j" %->% "i"))


dat %>%  filter(IsEdgePresent == TRUE, Stat == "Power") %>% 
  mutate(dist= factor(dist)) %>% 
  ggplot(data=.)+
  geom_line(aes(x=N,y=Stat_Val,group = dist, color = dist))+
  ggtitle(expression("Power -- Truth: j" %->% "i,   H_0: NOT j" %->% "i"))


#full and halfvariance
rbind(dat1 %>% mutate(type= "Half var"),dat %>% mutate(type="Norm var")) %>% 
  filter(N <= 100000) %>% 
  filter(IsEdgePresent == FALSE, Stat == "Power") %>% 
  mutate(dist= factor(dist)) %>% 
  ggplot(data=.)+
  geom_line(aes(x=N,y=Stat_Val,group = interaction(dist,type), color = dist,linetype =type))+
  ggtitle(expression("Power -- Truth: NOT j" %->% "i,   H_0: j" %->% "i"))




# 8 Nodes Two Confidence methods: 

SummaryDat <- readRDS(file="Data/Data_Testing_Multivariate_v2_2201130221.RDS")


dat <- SummaryDat %>%   
  unnest(cols=c("Test")) %>% 
  gather(c(Test_Bonferroni,Test_Simultaneous),key=Test_Type,value=Test) %>% 
  mutate(Stat = ifelse( (IsEdgePresent == "TRUE" & Hyp == 1) | (IsEdgePresent == FALSE & Hyp == -1 ),"Level","Power")) %>% 
  mutate(dist = case_when(dist < 0 ~ "<0", dist == 0 ~ "0", dist > 0 ~ ">0")) %>% 
  group_by(N,p,TreeType,IsEdgePresent,dist,Hyp,Stat,Test_Type) %>% 
  summarise( from=mean(from),to=mean(to),Stat_Val = mean(Test), n = n()) %>%  
  arrange(IsEdgePresent,N,Stat,dist) %>% 
  print(n=200)

dat %>%  filter(IsEdgePresent == FALSE, Stat == "Power") %>% 
  mutate(dist= factor(dist)) %>% 
  ggplot(data=.)+
  geom_line(aes(x=N,y=Stat_Val,group = interaction(dist,Test_Type), color = dist,linetype=Test_Type))+
  ggtitle(expression("Power -- Truth: NOT j" %->% "i,   H_0: j" %->% "i"))


dat %>%  filter(IsEdgePresent == TRUE, Stat == "Power") %>% 
  mutate(dist= factor(dist)) %>% 
  ggplot(data=.)+
  geom_line(aes(x=N,y=Stat_Val,group = dist, color = dist))+
  ggtitle(expression("Power -- Truth: j" %->% "i,   H_0: NOT j" %->% "i"))



# 8 nodes 50k samples max Two Confidence Two Test Methods:



SummaryDat <- readRDS(file="Data/Data_Testing_Multivariate_v2_2201131523.RDS")


dat <- SummaryDat %>%   
  unnest(cols=c("Test")) %>% 
  gather(c(Test_Bonferroni,Test_Simultaneous,Test2_Bonferroni,Test2_Simultaneous),key=Test_Type,value=Test) %>% 
  mutate(Stat = ifelse( (IsEdgePresent == "TRUE" & Hyp == 1) | (IsEdgePresent == FALSE & Hyp == -1 ),"Level","Power")) %>% 
  mutate(dist = case_when(dist < 0 ~ "<0", dist == 0 ~ "0", dist > 0 ~ ">0")) %>% 
  group_by(N,p,TreeType,IsEdgePresent,dist,Hyp,Stat,Test_Type) %>% 
  summarise( from=mean(from),to=mean(to),Stat_Val = mean(Test), n = n()) %>%  
  arrange(IsEdgePresent,N,Stat,dist) %>% 
  print(n=200)

p1 <- dat %>%  filter(IsEdgePresent == FALSE, Stat == "Power") %>% 
  mutate(temp = stringr::str_split(Test_Type, '_')) %>% 
  rowwise() %>% 
  mutate(TestMethod= ifelse(unlist(temp)[1] == "Test","Old","New"), ConfidenceMethod = unlist(temp)[2]) %>% 
  ungroup() %>% 
  mutate(dist= as.factor(dist)) %>% 
  ggplot(data=.)+
  geom_line(aes(x=N,y=Stat_Val,group = interaction(dist,TestMethod,ConfidenceMethod), color = dist,linetype=ConfidenceMethod,size=TestMethod))+
  ggtitle(expression("Power -- Truth: NOT j" %->% "i,   H_0: j" %->% "i"))+
  ylab("Power")+
  theme(legend.position="none")


p2 <- dat %>%  filter(IsEdgePresent == TRUE, Stat == "Power") %>% 
  mutate(temp = stringr::str_split(Test_Type, '_')) %>% 
  rowwise() %>% 
  mutate(TestMethod= ifelse(unlist(temp)[1] == "Test","Old","New"), ConfidenceMethod = unlist(temp)[2]) %>% 
  ungroup() %>% 
  mutate(dist= as.factor(dist)) %>% 
  ggplot(data=.)+
  geom_line(aes(x=N,y=Stat_Val,group = interaction(dist,TestMethod,ConfidenceMethod), color = dist,linetype=ConfidenceMethod,size=TestMethod))+
  ggtitle(expression("Power -- Truth: j" %->% "i,   H_0: NOT j" %->% "i"))+
  ylab("Power")



p3 <- dat %>%  filter(IsEdgePresent == FALSE, Stat == "Level") %>% 
  mutate(temp = stringr::str_split(Test_Type, '_')) %>% 
  rowwise() %>% 
  mutate(TestMethod= ifelse(unlist(temp)[1] == "Test","Old","New"), ConfidenceMethod = unlist(temp)[2]) %>% 
  ungroup() %>% 
  mutate(dist= as.factor(dist)) %>% 
  ggplot(data=.)+
  geom_line(aes(x=N,y=Stat_Val,group = interaction(dist,TestMethod,ConfidenceMethod), color = dist,linetype=ConfidenceMethod,size=TestMethod))+
  ggtitle(expression("Level -- Truth: NOT j" %->% "i,   H_0: NOT j" %->% "i")) +
  ylab("Level")+
  theme(legend.position="none")


p4 <- dat %>%  filter(IsEdgePresent == TRUE, Stat == "Level") %>% 
  mutate(temp = stringr::str_split(Test_Type, '_')) %>% 
  rowwise() %>% 
  mutate(TestMethod= ifelse(unlist(temp)[1] == "Test","Old","New"), ConfidenceMethod = unlist(temp)[2]) %>% 
  ungroup() %>% 
  mutate(dist= as.factor(dist)) %>% 
  ggplot(data=.)+
  geom_line(aes(x=N,y=Stat_Val,group = interaction(dist,TestMethod,ConfidenceMethod), color = dist,linetype=ConfidenceMethod,size=TestMethod))+
  ggtitle(expression("Level -- Truth: j" %->% "i,   H_0:  j" %->% "i"))+
  ylab("Level")

ggarrange(p1,p2,p3,p4,ncol=2, nrow=2, common.legend = TRUE, legend="right")


# 2 and 4 nodes 20k samples max Two Confidence Two Test Methods:  

SummaryDat <- readRDS(file="Data/Data_Testing_Multivariate_v2_2201141350.RDS")


dat <- SummaryDat %>%   
  unnest(cols=c("Test")) %>% 
  gather(c(Test_Bonferroni,Test2_Bonferroni),key=Test_Type,value=Test) %>% 
  mutate(Stat = ifelse( (IsEdgePresent == "TRUE" & Hyp == 1) | (IsEdgePresent == FALSE & Hyp == -1 ),"Level","Power")) %>% 
  mutate(dist = case_when(dist < 0 ~ "<0", dist == 0 ~ "0", dist > 0 ~ ">0")) %>% 
  group_by(N,p,TreeType,IsEdgePresent,dist,Hyp,Stat,Test_Type) %>% 
  summarise( from=mean(from),to=mean(to),Stat_Val = mean(Test), n = n()) %>%  
  arrange(IsEdgePresent,N,Stat,dist) %>% 
  mutate(p =factor(p)) %>% 
  print(n=200)

p1 <- dat %>%  filter(IsEdgePresent == FALSE, Stat == "Power") %>% 
  mutate(temp = stringr::str_split(Test_Type, '_')) %>% 
  rowwise() %>% 
  mutate(TestMethod= ifelse(unlist(temp)[1] == "Test","Old","New"), ConfidenceMethod = unlist(temp)[2]) %>% 
  ungroup() %>% 
  mutate(dist= as.factor(dist)) %>% 
  ggplot(data=.)+
  geom_line(aes(x=N,y=Stat_Val,group = interaction(dist,TestMethod,p), color = dist,linetype=p,size=TestMethod))+
  ggtitle(expression("Power -- Truth: NOT j" %->% "i,   H_0: j" %->% "i"))+
  ylab("Power")+
  theme(legend.position="none")


p2 <- dat %>%  filter(IsEdgePresent == TRUE, Stat == "Power") %>% 
  mutate(temp = stringr::str_split(Test_Type, '_')) %>% 
  rowwise() %>% 
  mutate(TestMethod= ifelse(unlist(temp)[1] == "Test","Old","New"), ConfidenceMethod = unlist(temp)[2]) %>% 
  ungroup() %>% 
  mutate(dist= as.factor(dist)) %>% 
  ggplot(data=.)+
  geom_line(aes(x=N,y=Stat_Val,group = interaction(dist,TestMethod,p), color = dist,linetype=p,size=TestMethod))+
  ggtitle(expression("Power -- Truth: j" %->% "i,   H_0: NOT j" %->% "i"))+
  ylab("Power")



p3 <- dat %>%  filter(IsEdgePresent == FALSE, Stat == "Level") %>% 
  mutate(temp = stringr::str_split(Test_Type, '_')) %>% 
  rowwise() %>% 
  mutate(TestMethod= ifelse(unlist(temp)[1] == "Test","Old","New"), ConfidenceMethod = unlist(temp)[2]) %>% 
  ungroup() %>% 
  mutate(dist= as.factor(dist)) %>% 
  ggplot(data=.)+
  geom_line(aes(x=N,y=Stat_Val,group = interaction(dist,TestMethod,p), color = dist,linetype=p,size=TestMethod))+
  ggtitle(expression("Level -- Truth: NOT j" %->% "i,   H_0: NOT j" %->% "i")) +
  ylab("Level")+
  theme(legend.position="none")


p4 <- dat %>%  filter(IsEdgePresent == TRUE, Stat == "Level") %>% 
  mutate(temp = stringr::str_split(Test_Type, '_')) %>% 
  rowwise() %>% 
  mutate(TestMethod= ifelse(unlist(temp)[1] == "Test","Old","New"), ConfidenceMethod = unlist(temp)[2]) %>% 
  ungroup() %>% 
  mutate(dist= as.factor(dist)) %>% 
  ggplot(data=.)+
  geom_line(aes(x=N,y=Stat_Val,group = interaction(dist,TestMethod,p), color = dist,linetype=p,size=TestMethod))+
  ggtitle(expression("Level -- Truth: j" %->% "i,   H_0:  j" %->% "i"))+
  ylab("Level")

ggarrange(p1,p2,p3,p4,ncol=2, nrow=2, common.legend = TRUE, legend="right")


# 2 Nodes
SummaryDat <- readRDS(file="Data/Data_Bivariate_Testing_2201071517.RDS")

SummaryDat %>% unnest(cols=c("Test")) %>% 
  mutate(Reject = ifelse(Test==TRUE, 0, 1)) %>% 
  group_by(N,p,TreeType,HypothesisType,HypothesisLocation) %>% 
  summarise( from=mean(from),to=mean(to),Stat_Val = mean(Reject)) %>%  
  mutate(Stat = ifelse(HypothesisType == "Missing","Power","Level")) %>% 
  arrange(HypothesisType,HypothesisLocation) %>%  
  select(N,p,TreeType,HypothesisType,HypothesisLocation,from,to,Stat,Stat_Val) %>% 
  print(n=100)





################# TABLE ###################

# 2 and 4 nodes 20k samples max Two Confidence Two Test Methods:  

SummaryDat_2_4 <- readRDS(file="Data/Data_Testing_Multivariate_v2_2201141350.RDS")


stat_2_4_total <- SummaryDat_2_4 %>% 
  unnest(cols=c("Test")) %>% 
  gather(c(Test_Bonferroni,Test2_Bonferroni),key=Test_Type,value=Test) %>% 
  filter(Test_Type == "Test2_Bonferroni")  %>% 
  select(-Test_Type) %>% 
  mutate(Stat = ifelse( (IsEdgePresent == "TRUE" & Hyp == 1) | (IsEdgePresent == FALSE & Hyp == -1 ),"Level","Power")) %>% 
  #filter(Stat == "Power") %>% 
  mutate(dist = case_when(dist < 0 ~ "<0", dist == 0 ~ "0", dist > 0 ~ ">0")) %>%
  group_by(.,N,p,TreeType,IsEdgePresent,Hyp,Stat) %>% 
  summarise(Stat_Val = mean(Test)) %>% 
  mutate(col = paste(Stat,IsEdgePresent,sep="_")) %>% 
  ungroup %>% 
  select(-IsEdgePresent,-Stat,-Hyp) %>% 
  pivot_wider(names_from = col, values_from = Stat_Val) %>% 
  arrange(p,N)

Power_2_4_bydist <- SummaryDat_2_4 %>% 
  unnest(cols=c("Test")) %>% 
  gather(c(Test_Bonferroni,Test2_Bonferroni),key=Test_Type,value=Test) %>% 
  filter(Test_Type == "Test2_Bonferroni")  %>% 
  select(-Test_Type) %>% 
  mutate(Stat = ifelse( (IsEdgePresent == "TRUE" & Hyp == 1) | (IsEdgePresent == FALSE & Hyp == -1 ),"Level","Power")) %>% 
  filter(Stat == "Power") %>% 
  mutate(dist = case_when(dist < 0 ~ "<0", dist == 0 ~ "0", dist > 0 ~ ">0")) %>%
  group_by(.,N,p,TreeType,IsEdgePresent,dist,Hyp,Stat) %>% 
  summarise(Stat_Val = mean(Test)) %>% 
  mutate(col = paste(Stat,IsEdgePresent,dist,sep="_")) %>% 
  ungroup %>% 
  select(-IsEdgePresent,-Stat,-Hyp,-dist) %>% 
  pivot_wider(names_from = col, values_from = Stat_Val) %>% 
  arrange(p,N)

stat_2_4_tab <- left_join(stat_2_4_total,Power_2_4_bydist) %>% 
  select(p,N,'Power_FALSE_<0','Power_FALSE_>0','Power_FALSE_0','Power_FALSE','Power_TRUE','Level_TRUE','Level_FALSE')


# 8 nodes
SummaryDat_8 <- readRDS(file="Data/Data_Testing_Multivariate_v2_2201142102.RDS")



stat_8_total <- SummaryDat_8 %>% 
  unnest(cols=c("Test")) %>% 
  gather(c(Test2_Bonferroni),key=Test_Type,value=Test) %>% 
  filter(Test_Type == "Test2_Bonferroni")  %>% 
  select(-Test_Type) %>% 
  mutate(Stat = ifelse( (IsEdgePresent == "TRUE" & Hyp == 1) | (IsEdgePresent == FALSE & Hyp == -1 ),"Level","Power")) %>% 
  #filter(Stat == "Power") %>% 
  mutate(dist = case_when(dist < 0 ~ "<0", dist == 0 ~ "0", dist > 0 ~ ">0")) %>%
  group_by(.,N,p,TreeType,IsEdgePresent,Hyp,Stat) %>% 
  summarise(Stat_Val = mean(Test)) %>% 
  mutate(col = paste(Stat,IsEdgePresent,sep="_")) %>% 
  ungroup %>% 
  select(-IsEdgePresent,-Stat,-Hyp) %>% 
  pivot_wider(names_from = col, values_from = Stat_Val) %>% 
  arrange(p,N)

Power_8_bydist <- SummaryDat_8 %>% 
  unnest(cols=c("Test")) %>% 
  gather(c(Test2_Bonferroni),key=Test_Type,value=Test) %>% 
  filter(Test_Type == "Test2_Bonferroni")  %>% 
  select(-Test_Type) %>% 
  mutate(Stat = ifelse( (IsEdgePresent == "TRUE" & Hyp == 1) | (IsEdgePresent == FALSE & Hyp == -1 ),"Level","Power")) %>% 
  filter(Stat == "Power") %>% 
  mutate(dist = case_when(dist < 0 ~ "<0", dist == 0 ~ "0", dist > 0 ~ ">0")) %>%
  group_by(.,N,p,TreeType,IsEdgePresent,dist,Hyp,Stat) %>% 
  summarise(Stat_Val = mean(Test)) %>% 
  mutate(col = paste(Stat,IsEdgePresent,dist,sep="_")) %>% 
  ungroup %>% 
  select(-IsEdgePresent,-Stat,-Hyp,-dist) %>% 
  pivot_wider(names_from = col, values_from = Stat_Val) %>% 
  arrange(p,N)

stat_8_tab <- left_join(stat_8_total,Power_8_bydist) %>% 
  select(p,N,'Power_FALSE_<0','Power_FALSE_>0','Power_FALSE_0','Power_FALSE','Power_TRUE','Level_TRUE','Level_FALSE')


#16 nodes

SummaryDat_16 <- readRDS(file="Data/Data_Testing_Multivariate_v2_2201141652.RDS")



stat_16_total <- SummaryDat_16 %>% 
  unnest(cols=c("Test")) %>% 
  gather(c(Test2_Bonferroni),key=Test_Type,value=Test) %>% 
  filter(Test_Type == "Test2_Bonferroni")  %>% 
  select(-Test_Type) %>% 
  mutate(Stat = ifelse( (IsEdgePresent == "TRUE" & Hyp == 1) | (IsEdgePresent == FALSE & Hyp == -1 ),"Level","Power")) %>% 
  #filter(Stat == "Power") %>% 
  mutate(dist = case_when(dist < 0 ~ "<0", dist == 0 ~ "0", dist > 0 ~ ">0")) %>%
  group_by(.,N,p,TreeType,IsEdgePresent,Hyp,Stat) %>% 
  summarise(Stat_Val = mean(Test)) %>% 
  mutate(col = paste(Stat,IsEdgePresent,sep="_")) %>% 
  ungroup %>% 
  select(-IsEdgePresent,-Stat,-Hyp) %>% 
  pivot_wider(names_from = col, values_from = Stat_Val) %>% 
  arrange(p,N)

Power_16_bydist <- SummaryDat_16 %>% 
  unnest(cols=c("Test")) %>% 
  gather(c(Test2_Bonferroni),key=Test_Type,value=Test) %>% 
  filter(Test_Type == "Test2_Bonferroni")  %>% 
  select(-Test_Type) %>% 
  mutate(Stat = ifelse( (IsEdgePresent == "TRUE" & Hyp == 1) | (IsEdgePresent == FALSE & Hyp == -1 ),"Level","Power")) %>% 
  filter(Stat == "Power") %>% 
  mutate(dist = case_when(dist < 0 ~ "<0", dist == 0 ~ "0", dist > 0 ~ ">0")) %>%
  group_by(.,N,p,TreeType,IsEdgePresent,dist,Hyp,Stat) %>% 
  summarise(Stat_Val = mean(Test)) %>% 
  mutate(col = paste(Stat,IsEdgePresent,dist,sep="_")) %>% 
  ungroup %>% 
  select(-IsEdgePresent,-Stat,-Hyp,-dist) %>% 
  pivot_wider(names_from = col, values_from = Stat_Val) %>% 
  arrange(p,N)

stat_16_tab <- left_join(stat_16_total,Power_16_bydist) %>% 
  select(p,N,'Power_FALSE_<0','Power_FALSE_>0','Power_FALSE_0','Power_FALSE','Power_TRUE','Level_TRUE','Level_FALSE')


#Samlet  200 rep

SummaryDat <- readRDS(file="Data/Data_Testing_Multivariate_v2_2201142322.RDS")


