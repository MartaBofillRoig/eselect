##################################################################################
# Research project Vienna - Bcn
# Endpoint selection and sample size reassessment for composite binary endpoints
# Results Simulation study 
##################################################################################

rm(list = ls())

setwd("H:/Code_PAPER_biost")
# setwd("~/Code_PAPER_biost")

library(ggplot2)
library(gridExtra)
library(ggpubr)
library(tidyverse) 

##################################################################################
# Results under H0 - Study on the significance level
##################################################################################

# Under H0
load("H:/Code_PAPER_biost/results/results_H0True.RData") 
# load("~/Code_PAPER_biost/results/results_H0True.RData") 
dataset_H0True = dataset 

names(dataset_H0True)[names(dataset_H0True)=="Test_Power_ES"] <- 'Test_Reject_ES'
names(dataset_H0True)[names(dataset_H0True)=="Test_Power_ES_SS"] <- 'Test_Reject_ES_SS'
names(dataset_H0True)[names(dataset_H0True)=="Test_Power_ES_ub"] <- 'Test_Reject_ES_ub'
names(dataset_H0True)[names(dataset_H0True)=="Test_Power_ES_ubSS"] <- 'Test_Reject_ES_ubSS'
names(dataset_H0True)[names(dataset_H0True)=="Test_Power_CE"] <- 'Test_Reject_CE'
names(dataset_H0True)[names(dataset_H0True)=="Test_Power_RE"] <- 'Test_Reject_RE'

signlevel_data <- data.frame(signlevel=c(dataset_H0True$Test_Reject_ES_SS,
                                         dataset_H0True$Test_Reject_ES_ubSS,
                                         dataset_H0True$Test_Reject_CE,dataset_H0True$Test_Reject_RE),
                             Endpoint=c(rep("AD",length(dataset_H0True$Test_Reject_ES_SS)),
                                        rep("AD (unblinded)",length(dataset_H0True$Test_Reject_ES_ubSS)),
                                        rep("CD",length(dataset_H0True$Test_Reject_CE)),
                                        rep("RD",length(dataset_H0True$Test_Reject_RE)))
)

ggplot(signlevel_data, aes(x=Endpoint, y=signlevel)) + geom_boxplot() + labs(y = "Significance level", x="Endpoint")
ggsave("H0True_unblinded.pdf",width = 100, height = 100, units = "mm")

subset_h0 <- data.frame(Test_Reject_ES_SS=dataset_H0True$Test_Reject_ES_SS,
                        Test_Reject_ES_ubSS=dataset_H0True$Test_Reject_ES_ubSS,
                        Test_Reject_CE=dataset_H0True$Test_Reject_CE,
                        Test_Reject_RE=dataset_H0True$Test_Reject_RE
)
summary(subset_h0)


summary(subset(dataset_H0True,dataset_H0True$p_init==1))
summary(subset(dataset_H0True,dataset_H0True$p_init==0.5)) 

##################################################################################
# Results under H1 - Study on the power
# with sample size reassessment
##################################################################################
rm(list = ls())

# Under H1
load("H:/Code_PAPER_biost/results/results_H0False_withss.RData")
# load("~/Code_PAPER_biost/results/results_H0False_withss.RData")
dataset_H0False_ss = dataset 
# rm(dataset)

# selecting only those scenarios with interim analysis
dataset_H0False_ss = subset(dataset,dataset$p_init==0.5) 
rm(dataset)

# 
summary(dataset_H0False_ss)

# 
power_data <- data.frame(Power=c(dataset_H0False_ss$Test_Power_ES_SS,dataset_H0False_ss$Test_Power_ES_ubSS,
                                 dataset_H0False_ss$Test_Power_CE,dataset_H0False_ss$Test_Power_RE),
                         Endpoint=c(rep("AD",length(dataset_H0False_ss$Test_Power_ES_SS)),
                                    rep("AD (unblinded)",length(dataset_H0False_ss$Test_Power_ES_ubSS)),
                                    rep("CD",length(dataset_H0False_ss$Test_Power_CE)),
                                    rep("RD",length(dataset_H0False_ss$Test_Power_RE)))
)
ggplot(power_data, aes(x=Endpoint, y=Power)) + geom_boxplot()
# ggsave("power_H0False_withssunblinded.pdf",width = 100, height = 100, units = "mm")

# plot without unblinded 
power_data_blinded <- data.frame(Power=c(dataset_H0False_ss$Test_Power_ES_SS,
                                 # dataset_H0False_ss$Test_Power_ES_ubSS,
                                 dataset_H0False_ss$Test_Power_CE,dataset_H0False_ss$Test_Power_RE),
                         Design=c(rep("AD",length(dataset_H0False_ss$Test_Power_ES_SS)),
                                    # rep("SE_ub",length(dataset_H0False_ss$Test_Power_ES_ubSS)),
                                    rep("CD",length(dataset_H0False_ss$Test_Power_CE)),
                                    rep("RD",length(dataset_H0False_ss$Test_Power_RE)))
)
ggplot(power_data_blinded, aes(x=Design, y=Power)) + geom_boxplot()
# ggsave("power_H0False_withss.pdf",width = 100, height = 100, units = "mm")

################################################################################## 
p <- list() 
enum = 1 
it <- 1 

for(i in 1:max(dataset_H0False_ss$scenario)){
  sub=subset(dataset_H0False_ss,dataset_H0False_ss$scenario==i & dataset_H0False_ss$p_init==0.5) 
  univ=rbind(sub,sub)
  univ$testind=c(rep("RE",dim(sub)[1]),rep("CE",dim(sub)[1]))
  univ$test=ifelse(univ$testind=="RE",univ$Test_Power_RE,univ$Test_Power_CE)
  
  # 
  p[[enum]] <-ggplot(sub, aes(x=corr, y=Test_Power_ES_SS, color=as.factor(ss_decision)))  +
    geom_point(size=2)+ ggtitle(paste("Scenario", dataset_H0False_ss$scenario[it], "\n (p1,p2,OR1,OR2) \n=(", dataset_H0False_ss$p0_e1[it],",",dataset_H0False_ss$p0_e2[it],",",dataset_H0False_ss$OR1[it],",",dataset_H0False_ss$OR2[it],")"))+geom_point(size=2)  + labs(y = "Power (AD)", x="Correlation", color="Decision SS") + coord_cartesian(ylim = c(0.40, 1))+ geom_path()+ theme(plot.title = element_text(size=9),legend.position="bottom", legend.title = element_text(size = 6), legend.text = element_text(size = 6))
  # + theme(plot.title = element_text(size=9),legend.position = c(0.8, 0.2))
  
  p[[enum+1]] <-ggplot(sub, aes(x=corr, y=Test_Power_ES_ubSS, color=as.factor(ss_decision)))  +
    geom_point(size=2)+ ggtitle(paste("Scenario", dataset_H0False_ss$scenario[it], "\n (p1,p2,OR1,OR2) \n=(", dataset_H0False_ss$p0_e1[it],",",dataset_H0False_ss$p0_e2[it],",",dataset_H0False_ss$OR1[it],",",dataset_H0False_ss$OR2[it],")"))+geom_point(size=2)  + labs(y = "Power (AD, unblinded)", x="Correlation", color="Decision SS") + coord_cartesian(ylim = c(0.40, 1))+ geom_path()+ theme(plot.title = element_text(size=9),legend.position="bottom", legend.title = element_text(size = 6), legend.text = element_text(size = 6))
  # + theme(plot.title = element_text(size=9),legend.position = c(0.8, 0.2)) 
  
  #
  p[[enum+2]] <- ggplot(univ, aes(x=corr, y=test, color=as.factor(testind)))+
    geom_point(size=2)+ ggtitle(paste("Scenario", dataset_H0False_ss$scenario[it], "\n (p1,p2,OR1,OR2) \n=(", dataset_H0False_ss$p0_e1[it],",",dataset_H0False_ss$p0_e2[it],",",dataset_H0False_ss$OR1[it],",",dataset_H0False_ss$OR2[it],")"))+geom_point(size=2)  + labs(y = "Power (CD/RD)", x="Correlation", color="Decision") + coord_cartesian(ylim = c(0.40, 1))+ geom_path()+ theme(plot.title = element_text(size=9),legend.position="bottom",#legend.position = c(0.8, 0.2),
                  legend.title = element_text(size = 6), legend.text = element_text(size = 6))
  
  # stable <- data.frame(Corr=sub$corr, 
                       # SS=round(sub$ss_ratio,2),"%CE"=round(100*sub$decision_ES_SS,2),
                       # "%CE(u)"=round(100*sub$decision_ES_ubSS,2),check.names=FALSE)  
  stable <- data.frame("Correlation"=sub$corr, 
                       "Decision rule"=round(sub$ss_ratio,2),"% CE"=round(100*sub$decision_ES_SS,2),check.names=FALSE)  
  
  p[[enum+3]] <- ggtexttable(stable, rows = NULL, theme = ttheme(base_style = "default", base_size = 9))
  
  enum=enum+4 
  it <- it + dim(subset(dataset_H0False_ss,dataset_H0False_ss$scenario==i))
}
plot <- marrangeGrob(p,ncol=3,nrow=4,top=NULL) 
# ggsave(file="results_ss.pdf", plot, width = 210, height = 297, units = "mm") 
ggsave(file="results_ss2.pdf", plot, width = 210, height = 297, units = "mm") 
# 
################################################################################## 
p <- list() 
enum = 1 
it <- 1 

i=1
dataset=dataset_H0False_ss
# for(i in 1:max(dataset$scenario)){
sub=subset(dataset,dataset$scenario==i & dataset$p_init==0.5) 
univ=rbind(sub,sub)
univ$testind=c(rep("RE",dim(sub)[1]),rep("CE",dim(sub)[1]))
univ$test=ifelse(univ$testind=="RE",univ$Test_Power_RE,univ$Test_Power_CE)

# CE/RE
p[[enum]] <- ggplot(univ, aes(x=corr, y=test, color=as.factor(testind)))+
  geom_point(size=2)+ ggtitle(paste("Scenario:", dataset$scenario[it], " (p1,p2,OR1,OR2) =(", dataset$p0_e1[it],",",dataset$p0_e2[it],",",dataset$OR1[it],",",dataset$OR2[it],")"))+geom_point(size=2)  + labs(y = "Power (CD/RD)", x="Correlation", color="Endpoint") + coord_cartesian(ylim = c(0.60, 1))+ geom_path()+ theme(plot.title = element_text(size=10),legend.position="bottom", legend.title = element_text(size = 9), legend.text = element_text(size = 6))

# ES
p[[enum+1]] <-ggplot(sub, aes(x=corr, y=Test_Power_ES_SS, color=as.factor(ss_decision)))  +
  geom_point(size=2)+ ggtitle(paste("Scenario:", dataset$scenario[it], " (p1,p2,OR1,OR2) =(", dataset$p0_e1[it],",",dataset$p0_e2[it],",",dataset$OR1[it],",",dataset$OR2[it],")"))+geom_point(size=2)  + labs(y = "Power (AD)", x="Correlation", color="Endpoint") + coord_cartesian(ylim = c(0.60, 1))+ geom_path()+ theme(plot.title = element_text(size=10),legend.position="bottom", legend.title = element_text(size = 9), legend.text = element_text(size = 6))
# + theme(plot.title = element_text(size=9),legend.position = c(0.8, 0.2)) 

stable <- data.frame("Correlation"=sub$corr, 
                     "Decision rule"=round(sub$ss_ratio,2),"% CE"=round(100*sub$decision_ES_SS,2),check.names=FALSE)  

p[[enum+2]] <- ggtexttable(stable, rows = NULL, theme = ttheme(base_style = "default", base_size = 12))

enum=enum+3
it <- it + dim(subset(dataset,dataset$scenario==i))
# }
plot <- marrangeGrob(p,ncol=3,nrow=1,top=NULL) 
# ggsave(file="whatever.pdf", plot, width = 210, height = 297, units = "mm") 
ggsave(file="power_H0False_withss_example.pdf", plot, width = 350, height = 100, units = "mm") 

# windows(height = 4, width = 10)
# plot


##################################################################################
# Results under H1 - Study on the power
# without sample size reassessment
##################################################################################
rm(list = ls())

# Under H1
load("H:/Code_PAPER_biost/results/results_H0False_withoutss2.RData")
# load("~/Code_PAPER_biost/results/results_H0False_withoutss2.RData")

summary(subset(dataset,dataset$p_init==0.5)) 
summary(subset(dataset,dataset$p_init==1)) 

# selecting only those scenarios without interim analysis
dataset_H0False_ss = subset(dataset,dataset$p_init==1) 
# rm(dataset)

# 
summary(dataset_H0False_ss)

# 
power_data <- data.frame(Power=c(dataset_H0False_ss$Test_Power_ES_SS,dataset_H0False_ss$Test_Power_ES_ubSS,
                                 dataset_H0False_ss$Test_Power_CE,dataset_H0False_ss$Test_Power_RE),
                         Endpoint=c(rep("AD",length(dataset_H0False_ss$Test_Power_ES_SS)),
                                    rep("AD (unblinded)",length(dataset_H0False_ss$Test_Power_ES_ubSS)),
                                    rep("CD",length(dataset_H0False_ss$Test_Power_CE)),
                                    rep("RD",length(dataset_H0False_ss$Test_Power_RE)))
)

ggplot(power_data, aes(x=Endpoint, y=Power)) + geom_boxplot()
# ggsave("power_H0False_withoutssunblinded.pdf",width = 100, height = 100, units = "mm")

# plot without unblinded 
power_data_blinded <- data.frame(Power=c(dataset_H0False_ss$Test_Power_ES_SS,
                                         # dataset_H0False_ss$Test_Power_ES_ubSS,
                                         dataset_H0False_ss$Test_Power_CE,dataset_H0False_ss$Test_Power_RE),
                                 Design=c(rep("AD",length(dataset_H0False_ss$Test_Power_ES_SS)),
                                            # rep("SE_ub",length(dataset_H0False_ss$Test_Power_ES_ubSS)),
                                            rep("CD",length(dataset_H0False_ss$Test_Power_CE)),
                                            rep("RD",length(dataset_H0False_ss$Test_Power_RE)))
)
ggplot(power_data_blinded, aes(x=Design, y=Power)) + geom_boxplot()
# ggsave("power_H0False_withoutss.pdf",width = 100, height = 100, units = "mm")

################################################################################## 
p <- list() 
enum = 1 
it <- 1 

for(i in 1:max(dataset_H0False_ss$scenario)){
  sub=subset(dataset_H0False_ss,dataset_H0False_ss$scenario==i & dataset_H0False_ss$p_init==1.0) 
  univ=rbind(sub,sub)
  univ$testind=c(rep("RE",dim(sub)[1]),rep("CE",dim(sub)[1]))
  univ$test=ifelse(univ$testind=="RE",univ$Test_Power_RE,univ$Test_Power_CE)
  
  # 
  p[[enum]] <-ggplot(sub, aes(x=corr, y=Test_Power_ES_SS, color=as.factor(ss_decision)))  +
    geom_point(size=2)+ ggtitle(paste("Scenario", dataset_H0False_ss$scenario[it], "\n (p1,p2,OR1,OR2) \n=(", dataset_H0False_ss$p0_e1[it],",",dataset_H0False_ss$p0_e2[it],",",dataset_H0False_ss$OR1[it],",",dataset_H0False_ss$OR2[it],")"))+geom_point(size=2)  + labs(y = "Power (AD)", x="Correlation", color="Decision SS") + coord_cartesian(ylim = c(0.40, 1))+ geom_path()+ theme(plot.title = element_text(size=9),legend.position="bottom", legend.title = element_text(size = 6), legend.text = element_text(size = 6))
  # + theme(plot.title = element_text(size=9),legend.position = c(0.8, 0.2))
  
  p[[enum+1]] <-ggplot(sub, aes(x=corr, y=Test_Power_ES_ubSS, color=as.factor(ss_decision)))  +
    geom_point(size=2)+ ggtitle(paste("Scenario", dataset_H0False_ss$scenario[it], "\n (p1,p2,OR1,OR2) \n=(", dataset_H0False_ss$p0_e1[it],",",dataset_H0False_ss$p0_e2[it],",",dataset_H0False_ss$OR1[it],",",dataset_H0False_ss$OR2[it],")"))+geom_point(size=2)  + labs(y = "Power (AD, unblinded)", x="Correlation", color="Decision SS") + coord_cartesian(ylim = c(0.40, 1))+ geom_path()+ theme(plot.title = element_text(size=9),legend.position="bottom", legend.title = element_text(size = 6), legend.text = element_text(size = 6))
  # + theme(plot.title = element_text(size=9),legend.position = c(0.8, 0.2)) 
  
  #
  p[[enum+2]] <- ggplot(univ, aes(x=corr, y=test, color=as.factor(testind)))+
    geom_point(size=2)+ ggtitle(paste("Scenario", dataset_H0False_ss$scenario[it], "\n (p1,p2,OR1,OR2) \n=(", dataset_H0False_ss$p0_e1[it],",",dataset_H0False_ss$p0_e2[it],",",dataset_H0False_ss$OR1[it],",",dataset_H0False_ss$OR2[it],")"))+geom_point(size=2)  + labs(y = "Power (CD/RD)", x="Correlation", color="Decision") + coord_cartesian(ylim = c(0.40, 1))+ geom_path()+ theme(plot.title = element_text(size=9),legend.position="bottom",#legend.position = c(0.8, 0.2),
                                                                                                                                                                                                                                                                                                                                                                                                      legend.title = element_text(size = 6), legend.text = element_text(size = 6))
  
  
  stable <- data.frame("Correlation"=sub$corr, 
                       "Decision rule"=round(sub$ss_ratio,2),"% CE"=round(100*sub$decision_ES_SS,2),check.names=FALSE)  
  
  # stable <- data.frame(Corr=sub$corr, 
  #                      SS=round(sub$ss_ratio,2),"%CE"=round(100*sub$decision_ES_SS,2),
  #                      "%CE(u)"=round(100*sub$decision_ES_ubSS,2),check.names=FALSE)  
  
  p[[enum+3]] <- ggtexttable(stable, rows = NULL, theme = ttheme(base_style = "default", base_size = 9))
  
  enum=enum+4 
  it <- it + dim(subset(dataset_H0False_ss,dataset_H0False_ss$scenario==i))
}
plot <- marrangeGrob(p,ncol=3,nrow=4,top=NULL) 
# ggsave(file="results_withoutss.pdf", plot, width = 210, height = 297, units = "mm") 
ggsave(file="results_withoutss2.pdf", plot, width = 210, height = 297, units = "mm") 

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

# CE/RE
p[[enum]] <- ggplot(univ, aes(x=corr, y=test, color=as.factor(testind)))+
  geom_point(size=2)+ ggtitle(paste("Scenario:", dataset$scenario[it], " (p1,p2,OR1,OR2) =(", dataset$p0_e1[it],",",dataset$p0_e2[it],",",dataset$OR1[it],",",dataset$OR2[it],")"))+geom_point(size=2)  + labs(y = "Power (CD/RD)", x="Correlation", color="Endpoint") + coord_cartesian(ylim = c(0.60, 1))+ geom_path()+ theme(plot.title = element_text(size=10),legend.position="bottom", legend.title = element_text(size = 9), legend.text = element_text(size = 6))

# ES
p[[enum+1]] <-ggplot(sub, aes(x=corr, y=Test_Power_ES_SS, color=as.factor(ss_decision)))  +
  geom_point(size=2)+ ggtitle(paste("Scenario:", dataset$scenario[it], " (p1,p2,OR1,OR2) =(", dataset$p0_e1[it],",",dataset$p0_e2[it],",",dataset$OR1[it],",",dataset$OR2[it],")"))+geom_point(size=2)  + labs(y = "Power (AD)", x="Correlation", color="Endpoint") + coord_cartesian(ylim = c(0.60, 1))+ geom_path()+ theme(plot.title = element_text(size=10),legend.position="bottom", legend.title = element_text(size = 9), legend.text = element_text(size = 6))
# + theme(plot.title = element_text(size=9),legend.position = c(0.8, 0.2)) 

stable <- data.frame("Correlation"=sub$corr, 
                     "Decision rule"=round(sub$ss_ratio,2),"% CE"=round(100*sub$decision_ES_SS,2),check.names=FALSE)  

p[[enum+2]] <- ggtexttable(stable, rows = NULL, theme = ttheme(base_style = "default", base_size = 12))

enum=enum+3
it <- it + dim(subset(dataset,dataset$scenario==i))
# }
plot <- marrangeGrob(p,ncol=3,nrow=1,top=NULL) 
# ggsave(file="whatever.pdf", plot, width = 210, height = 297, units = "mm") 
ggsave(file="power_H0False_withoutss_example.pdf", plot, width = 350, height = 100, units = "mm") 


##################################################################################
# Effect only on the relevant endpoint
##################################################################################
rm(list = ls())

# Under H1
load("H:/Code_PAPER_biost/results/results_H0True_e1.RData")
# load("~/Code_PAPER_biost/results/results_H0True_e1.RData") 
dataset_H0False_ss = dataset 
rm(dataset)

# 
summary(dataset_H0False_ss)

# 
power_data <- data.frame(Power=c(dataset_H0False_ss$Test_Power_ES_SS,dataset_H0False_ss$Test_Power_ES_ubSS,
                                 dataset_H0False_ss$Test_Power_CE,dataset_H0False_ss$Test_Power_RE),
                         Design=c(rep("AD",length(dataset_H0False_ss$Test_Power_ES_SS)),
                                    rep("AD (unblinded)",length(dataset_H0False_ss$Test_Power_ES_ubSS)),
                                    rep("CE",length(dataset_H0False_ss$Test_Power_CE)),
                                    rep("RE",length(dataset_H0False_ss$Test_Power_RE)))
)

ggplot(power_data, aes(x=Design, y=Power)) + geom_boxplot()
ggsave("power_H0True_e1.pdf",width = 100, height = 100, units = "mm")

################################################################################## 
p <- list() 
enum = 1 
it <- 1 

for(i in 1:max(dataset_H0False_ss$scenario)){
  sub=subset(dataset_H0False_ss,dataset_H0False_ss$scenario==i & dataset_H0False_ss$p_init==1.0) 
  univ=rbind(sub,sub)
  univ$testind=c(rep("RE",dim(sub)[1]),rep("CE",dim(sub)[1]))
  univ$test=ifelse(univ$testind=="RE",univ$Test_Power_RE,univ$Test_Power_CE)
  
  # 
  p[[enum]] <-ggplot(sub, aes(x=corr, y=Test_Power_ES_SS, color=as.factor(ss_decision)))  +
    geom_point(size=2)+ ggtitle(paste("Scenario", dataset_H0False_ss$scenario[it], "\n (p1,p2,OR1,OR2) \n=(", dataset_H0False_ss$p0_e1[it],",",dataset_H0False_ss$p0_e2[it],",",dataset_H0False_ss$OR1[it],",",dataset_H0False_ss$OR2[it],")"))+geom_point(size=2)  + labs(y = "Power (AD)", x="Correlation", color="Decision SS") + coord_cartesian(ylim = c(0.60, 1))+ geom_path()+ theme(plot.title = element_text(size=9),legend.position="bottom", legend.title = element_text(size = 6), legend.text = element_text(size = 6))
  # + theme(plot.title = element_text(size=9),legend.position = c(0.8, 0.2))
  
  p[[enum+1]] <-ggplot(sub, aes(x=corr, y=Test_Power_ES_ubSS, color=as.factor(ss_decision)))  +
    geom_point(size=2)+ ggtitle(paste("Scenario", dataset_H0False_ss$scenario[it], "\n (p1,p2,OR1,OR2) \n=(", dataset_H0False_ss$p0_e1[it],",",dataset_H0False_ss$p0_e2[it],",",dataset_H0False_ss$OR1[it],",",dataset_H0False_ss$OR2[it],")"))+geom_point(size=2)  + labs(y = "Power (AD, unblinded)", x="Correlation", color="Decision SS") + coord_cartesian(ylim = c(0.60, 1))+ geom_path()+ theme(plot.title = element_text(size=9),legend.position="bottom", legend.title = element_text(size = 6), legend.text = element_text(size = 6))
  # + theme(plot.title = element_text(size=9),legend.position = c(0.8, 0.2)) 
  
  #
  p[[enum+2]] <- ggplot(univ, aes(x=corr, y=test, color=as.factor(testind)))+
    geom_point(size=2)+ ggtitle(paste("Scenario", dataset_H0False_ss$scenario[it], "\n (p1,p2,OR1,OR2) \n=(", dataset_H0False_ss$p0_e1[it],",",dataset_H0False_ss$p0_e2[it],",",dataset_H0False_ss$OR1[it],",",dataset_H0False_ss$OR2[it],")"))+geom_point(size=2)  + labs(y = "Power (CD/RD)", x="Correlation", color="Decision") + coord_cartesian(ylim = c(0.60, 1))+ geom_path()+ theme(plot.title = element_text(size=9),legend.position="bottom", legend.title = element_text(size = 6), legend.text = element_text(size = 6))
  
  stable <- data.frame(Corr=sub$corr, 
                       SS=round(sub$ss_ratio,2),"%CE"=round(100*sub$decision_ES_SS,2),
                       "%CE(u)"=round(100*sub$decision_ES_ubSS,2),check.names=FALSE)  
  
  p[[enum+3]] <- ggtexttable(stable, rows = NULL, theme = ttheme(base_style = "default", base_size = 9))
  
  enum=enum+4 
  it <- it + dim(subset(dataset_H0False_ss,dataset_H0False_ss$scenario==i))
}
plot <- marrangeGrob(p,ncol=3,nrow=4,top=NULL) 
# ggsave(file="whatever2.pdf", plot, width = 210, height = 297, units = "mm") 

##################################################################################

##################################################################################
# Effect only on the relevant endpoint
##################################################################################
rm(list = ls())

# Under H1
load("H:/Code_PAPER_biost/results/results_H0True_e2.RData")
# load("~/Code_PAPER_biost/results/results_H0True_e1.RData") 
dataset_H0True_e2 = dataset 
rm(dataset)

# 
summary(dataset_H0True_e2)

# 
power_data <- data.frame(Power=c(dataset_H0True_e2$Test_Power_ES_SS,dataset_H0True_e2$Test_Power_ES_ubSS,
                                 dataset_H0True_e2$Test_Power_CE,dataset_H0True_e2$Test_Power_RE),
                         Design=c(rep("AD",length(dataset_H0True_e2$Test_Power_ES_SS)),
                                    rep("AD (unblinded)",length(dataset_H0True_e2$Test_Power_ES_ubSS)),
                                    rep("CD",length(dataset_H0True_e2$Test_Power_CE)),
                                    rep("RD",length(dataset_H0True_e2$Test_Power_RE)))
)

ggplot(power_data, aes(x=Design, y=Power)) + geom_boxplot()
ggsave("power_H0True_e2.pdf",width = 100, height = 100, units = "mm")
