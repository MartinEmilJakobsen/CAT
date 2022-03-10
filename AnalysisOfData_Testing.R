library(scales)
library(ggpubr)
library(xtable)
library(kableExtra)
library(latex2exp)

# 2,4,6,8,16 nodes both test methods, 400 rep
SummaryDat <- readRDS(file="Data/Data_Testing_Multivariate_2201170534.RDS")


#Plot 
PlotData <- SummaryDat %>% 
  unnest(cols=c("Test")) %>% 
  rename(Exact = Test_Bonferroni,Asymptotic = Test2_Bonferroni) %>% 
  gather(c(Exact,Asymptotic),key='Query Scheme',value=Test) %>% 
  mutate(Stat = ifelse( (IsEdgePresent == "TRUE" & Hyp == 1) | (IsEdgePresent == FALSE & Hyp == -1 ),"Level","Power")) %>% 
  filter(Stat == "Power") %>% 
  group_by(.,N,p,TreeType,IsEdgePresent,Hyp,Stat,`Query Scheme`) %>% 
  summarise(Stat_Val = mean(Test)) %>% 
  mutate(col = paste(Stat,IsEdgePresent,sep="_")) %>% 
  ungroup %>% 
  select(-IsEdgePresent,-Stat,-Hyp) %>% 
  mutate(p = factor(p),
         n = N)

PlotData %>%  print(n=200)

PowerPlot <- ggplot(data=PlotData)+
  geom_line(aes(x=n,y=Stat_Val,group = interaction(p,`Query Scheme`), color = p,linetype=`Query Scheme`))+
  ylab("Power")+
  facet_wrap(~col,scales = "free")+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )


ggsave(
  "Plots/Testing_power.pdf",
  plot = PowerPlot,
  device = cairo_pdf,
  path = NULL,
  scale = 1,
  width = 210,
  height = 80,
  units = c("mm"),
  dpi = 500,
  limitsize = TRUE
)

#Table

stat_total <- SummaryDat %>% 
  unnest(cols=c("Test")) %>% 
  gather(c(Test_Bonferroni,Test2_Bonferroni),key=Test_Type,value=Test) %>% 
  mutate(Stat = ifelse( (IsEdgePresent == "TRUE" & Hyp == 1) | (IsEdgePresent == FALSE & Hyp == -1 ),"Level","Power")) %>% 
  #filter(Stat == "Power") %>% 
  mutate(dist = case_when(dist < 0 ~ "<0", dist == 0 ~ "0", dist > 0 ~ ">0")) %>%
  group_by(.,N,p,TreeType,IsEdgePresent,Hyp,Stat,Test_Type) %>% 
  summarise(Stat_Val = mean(Test)) %>% 
  mutate(col = paste(Stat,IsEdgePresent,sep="_")) %>% 
  ungroup %>% 
  select(-IsEdgePresent,-Stat,-Hyp) %>% 
  pivot_wider(names_from = col, values_from = Stat_Val) %>% 
  arrange(Test_Type,p,N)

Power_bydist <- SummaryDat %>% 
  unnest(cols=c("Test")) %>% 
  gather(c(Test_Bonferroni,Test2_Bonferroni),key=Test_Type,value=Test) %>% 
  #filter(Test_Type == "Test2_Bonferroni")  %>% 
  #select(-Test_Type) %>% 
  mutate(Stat = ifelse( (IsEdgePresent == "TRUE" & Hyp == 1) | (IsEdgePresent == FALSE & Hyp == -1 ),"Level","Power")) %>% 
  filter(Stat == "Power") %>% 
  mutate(dist = case_when(dist < 0 ~ "<0", dist == 0 ~ "0", dist > 0 ~ ">0")) %>%
  group_by(.,N,p,TreeType,IsEdgePresent,dist,Hyp,Stat,Test_Type) %>% 
  summarise(Stat_Val = mean(Test)) %>% 
  mutate(col = paste(Stat,IsEdgePresent,dist,sep="_")) %>% 
  ungroup %>% 
  select(-IsEdgePresent,-Stat,-Hyp,-dist) %>% 
  pivot_wider(names_from = col, values_from = Stat_Val) %>% 
  arrange(Test_Type,p,N)

Table <- left_join(stat_total,Power_bydist) %>% 
  select(Test_Type,p,N,'Power_FALSE_<0','Power_FALSE_>0','Power_FALSE_0','Power_FALSE','Power_TRUE','Level_TRUE','Level_FALSE') %>% 
  filter(Test_Type == "Test2_Bonferroni") %>% 
  select(-Test_Type)

#PrintTableToLatex:
options(knitr.kable.NA = '---')

  Table %>%  
    mutate_at(.vars=c('Power_FALSE_<0','Power_FALSE_>0','Power_FALSE_0','Power_FALSE','Power_TRUE'),.funs=function(x){floor(x* 100) / 100}) %>% 
    mutate_at(.vars=c('Level_TRUE','Level_FALSE'),.funs=function(x){ceiling(x* 100) / 100}) %>%   
    kable(format = 'latex', booktabs = TRUE) %>% 
    add_header_above(header = c("0" = 2, "Power" = 5,"Level"=2))
