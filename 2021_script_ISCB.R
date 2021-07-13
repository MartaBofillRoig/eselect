##################################################################################
# Research project Vienna - Bcn
# Endpoint selection and sample size reassessment for composite binary endpoints
# Results - ISCB
##################################################################################


rm(list = ls())

library(ggplot2) 
library(gridExtra)
library(ggpubr)
library(tidyverse) 

# 

load("C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/CBE_selection/results/results_H0False_withss.RData") 

summary(dataset)

# dataset=dataset[dataset$p_init==0.5,]

power_data <- data.frame(Power=c(dataset$Test_Power_ES_SS,
                                 dataset$Test_Power_CE,dataset$Test_Power_RE),
                         Endpoint=c(rep("Selected Endpoint",length(dataset$Test_Power_ES_SS)),
                                    rep("Composite Endpoint",length(dataset$Test_Power_CE)),
                                    rep("Relevant Endpoint",length(dataset$Test_Power_RE)))
)

ggplot_ss <- ggplot(power_data, aes(x=Endpoint, y=Power)) + geom_boxplot()

windows()
ggplot_ss

# 
################################################################################## 
p <- list() 
enum = 1 
it <- 1 

i=1
# for(i in 1:max(dataset$scenario)){
  sub=subset(dataset,dataset$scenario==i & dataset$p_init==1.0) 
  univ=rbind(sub,sub)
  univ$testind=c(rep("RE",dim(sub)[1]),rep("CE",dim(sub)[1]))
  univ$test=ifelse(univ$testind=="RE",univ$Test_Power_RE,univ$Test_Power_CE)
  
  # 
  p[[enum]] <-ggplot(sub, aes(x=corr, y=Test_Power_ES_SS, color=as.factor(ss_decision)))  +
    geom_point(size=2)+ ggtitle(paste("Scenario", dataset$scenario[it], "\n (p1,p2,OR1,OR2) \n=(", dataset$p0_e1[it],",",dataset$p0_e2[it],",",dataset$OR1[it],",",dataset$OR2[it],")"))+geom_point(size=2)  + labs(y = "Empirical Power (ES)", x="Correlation", color="Endpoint") + coord_cartesian(ylim = c(0.60, 1))+ geom_path()+ theme(plot.title = element_text(size=9),legend.position="bottom", legend.title = element_text(size = 8), legend.text = element_text(size = 6))
  # + theme(plot.title = element_text(size=9),legend.position = c(0.8, 0.2)) 
  
  #
  p[[enum+1]] <- ggplot(univ, aes(x=corr, y=test, color=as.factor(univ$testind)))+
    geom_point(size=2)+ ggtitle(paste("Scenario", dataset$scenario[it], "\n (p1,p2,OR1,OR2) \n=(", dataset$p0_e1[it],",",dataset$p0_e2[it],",",dataset$OR1[it],",",dataset$OR2[it],")"))+geom_point(size=2)  + labs(y = "Empirical Power (CE/RE)", x="Correlation", color="Endpoint") + coord_cartesian(ylim = c(0.60, 1))+ geom_path()+ theme(plot.title = element_text(size=9),legend.position="bottom", legend.title = element_text(size = 8), legend.text = element_text(size = 6))
  
  stable <- data.frame("Correlation"=sub$corr, 
                       "Decision rule"=round(sub$ss_ratio,2),"% Composite Endpoint"=round(100*sub$decision_ES_SS,2),check.names=FALSE)  
  
  p[[enum+2]] <- ggtexttable(stable, rows = NULL, theme = ttheme(base_style = "default", base_size = 9))
  
  enum=enum+3
  it <- it + dim(subset(dataset,dataset$scenario==i))
# }
plot <- marrangeGrob(p,ncol=3,nrow=1,top=NULL) 
# ggsave(file="whatever.pdf", plot, width = 210, height = 297, units = "mm") 

windows()
plot

################################################################################## 

load("C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/CBE_selection/20210-12-03b-with ARE/results/results_H0False_withoutss_new.RData")

power_data <- data.frame(Power=c(dataset$Test_Power_ES_SS,
                                 dataset$Test_Power_CE,dataset$Test_Power_RE),
                         Endpoint=c(rep("Selected Endpoint",length(dataset$Test_Power_ES_SS)),
                                    rep("Composite Endpoint",length(dataset$Test_Power_CE)),
                                    rep("Relevant Endpoint",length(dataset$Test_Power_RE)))
)

ggplot_ss <- ggplot(power_data, aes(x=Endpoint, y=Power)) + geom_boxplot()

windows()
ggplot_ss

################################################################################## 
p <- list() 
enum = 1 
it <- 1 

i=1
# for(i in 1:max(dataset$scenario)){
sub=subset(dataset,dataset$scenario==i & dataset$p_init==1.0) 
univ=rbind(sub,sub)
univ$testind=c(rep("RE",dim(sub)[1]),rep("CE",dim(sub)[1]))
univ$test=ifelse(univ$testind=="RE",univ$Test_Power_RE,univ$Test_Power_CE)

# 
p[[enum]] <-ggplot(sub, aes(x=corr, y=Test_Power_ES_SS, color=as.factor(ss_decision)))  +
  geom_point(size=2)+ ggtitle(paste("Scenario", dataset$scenario[it], "\n (p1,p2,OR1,OR2) \n=(", dataset$p0_e1[it],",",dataset$p0_e2[it],",",dataset$OR1[it],",",dataset$OR2[it],")"))+geom_point(size=2)  + labs(y = "Empirical Power (ES)", x="Correlation", color="Endpoint") + coord_cartesian(ylim = c(0.60, 1))+ geom_path()+ theme(plot.title = element_text(size=9),legend.position="bottom", legend.title = element_text(size = 8), legend.text = element_text(size = 6))
# + theme(plot.title = element_text(size=9),legend.position = c(0.8, 0.2)) 

#
p[[enum+1]] <- ggplot(univ, aes(x=corr, y=test, color=as.factor(univ$testind)))+
  geom_point(size=2)+ ggtitle(paste("Scenario", dataset$scenario[it], "\n (p1,p2,OR1,OR2) \n=(", dataset$p0_e1[it],",",dataset$p0_e2[it],",",dataset$OR1[it],",",dataset$OR2[it],")"))+geom_point(size=2)  + labs(y = "Empirical Power (CE/RE)", x="Correlation", color="Endpoint") + coord_cartesian(ylim = c(0.60, 1))+ geom_path()+ theme(plot.title = element_text(size=9),legend.position="bottom", legend.title = element_text(size = 8), legend.text = element_text(size = 6))

stable <- data.frame("Correlation"=sub$corr, 
                     "Decision rule"=round(sub$ss_ratio,2),"% Composite Endpoint"=round(100*sub$decision_ES_SS,2),check.names=FALSE)  

p[[enum+2]] <- ggtexttable(stable, rows = NULL, theme = ttheme(base_style = "default", base_size = 9))

enum=enum+3
it <- it + dim(subset(dataset,dataset$scenario==i))
# }
plot <- marrangeGrob(p,ncol=3,nrow=1,top=NULL) 
# ggsave(file="whatever.pdf", plot, width = 210, height = 297, units = "mm") 

windows()
plot

