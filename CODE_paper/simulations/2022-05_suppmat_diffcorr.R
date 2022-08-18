##################################################################################
# Research project Vienna - Bcn
# Endpoint selection and sample size reassessment for composite binary endpoints
# Simulation study for supp material
# Evaluating the power of the proposed design under diff correlations
##################################################################################

#########################################
rm(list = ls())

setwd("C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/CBE_selection/CODE_paper/simulations")
# setwd("H:/Code_PAPER_biost")
source('C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/CBE_selection/CODE_paper/simulations/aux_functions.R') 
source('C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/CBE_selection/eselect/R/eselect.R')
# source('H:/Code_PAPER_biost/eselect/R/eselectsim.R') 
# source('H:/Code_PAPER_biost/eselect/R/eselect_ub.R')
# source('H:/Code_PAPER_biost/eselect/R/eselectsim_ub.R') 

#########################################
# Preamble
#########################################


# install.packages(c("CompAREdesign","tidyverse","tidyr","plyr","devtools","ggplot2","gridExtra","ggpubr"))
# install.packages("tidyverse")
# install.packages("tidyr")
# install.packages("plyr")
# install.packages("devtools")
# install.packages("CompAREdesign")
# install.packages("ggplot2")
# install.packages("gridExtra")
# install.packages("ggpubr")

library(tidyverse)
library(tidyr)
library(plyr)
library(devtools) 
library(CompAREdesign)
library(ggplot2)
library(gridExtra)
library(ggpubr) 


f_dcorr <- function(ss_arm,p0_e1,OR1,p0_e2,OR2,p0_ce,p1_ce,p_init=1,criteria="SS",alpha=0.05,beta=0.2){
  
  n_init=ss_arm
  ss_arm=round(n_init*p_init)
  total_ss = ss_arm*2 
  
  p1_e1 = (OR1*p0_e1/(1-p0_e1))/(1+(OR1*p0_e1/(1-p0_e1)))
  p1_e2 = (OR2*p0_e2/(1-p0_e2))/(1+(OR2*p0_e2/(1-p0_e2)))
  
  # control group
  sm0 = f_sim(samplesize=ss_arm,p_e1=p0_e1,p_e2=p0_e2,p_ce=p0_ce)
  
  # intervention group
  sm1 = f_sim(samplesize=ss_arm,p_e1=p1_e1,p_e2=p1_e2,p_ce=p1_ce)
  
  # pooled sample
  sm = sm0 + sm1
  
  selection = eselect(db=sm,p0_e1=p0_e1,OR1=OR1,p0_e2=p0_e2,OR2=OR2,criteria=criteria,alpha=alpha,beta=beta)
  
  if(selection$Decision == 1){
    
    ss_arm_old = ss_arm
    ss_arm = n_init
    
    # control group
    sm0_add = f_sim(samplesize=round(ss_arm - ss_arm_old),p_e1=p0_e1,p_e2=p0_e2,p_ce=p0_ce)
    sm0 = sm0 + sm0_add
    # intervention group
    sm1_add = f_sim(samplesize=round(ss_arm - ss_arm_old),p_e1=p1_e1,p_e2=p1_e2,p_ce=p1_ce)
    sm1 = sm1 + sm1_add
    
    phat_group1 = 1-(sm1[4])/ss_arm
    phat_group0 = 1-(sm0[4])/ss_arm
    
    # test odds ratio CE
    TestOR_unpooled = test_f(OR=(phat_group1/(1-phat_group1))/(phat_group0/(1-phat_group0)),p0=phat_group0,n=ss_arm)
    
  }else{
    
    ss_arm_old = ss_arm
    ss_arm = n_init
    
    # control group
    sm0_add = f_sim(samplesize=round(ss_arm - ss_arm_old),p_e1=p0_e1,p_e2=p0_e2,p_ce=p0_ce)
    sm0 = sm0 + sm0_add
    # intervention group
    sm1_add = f_sim(samplesize=round(ss_arm - ss_arm_old),p_e1=p1_e1,p_e2=p1_e2,p_ce=p1_ce)
    sm1 = sm1 + sm1_add
    
    phat_group1 = (sm1[1]+sm1[2])/ss_arm
    phat_group0 = (sm0[1]+sm0[2])/ss_arm
    
    # test odds ratio RE
    TestOR_unpooled = test_f(OR=(phat_group1/(1-phat_group1))/(phat_group0/(1-phat_group0)),p0=phat_group0,n=ss_arm)
  }
  
  output = list(Test=TestOR_unpooled, Decision=selection$Decision)
  
  return(output)
}

#########################################
# Define the set of scenarios
#########################################

# Scenarios
p0_e1 = c(0.1)
p0_e2 = c(0.25)
OR1 = c(0.6) #old_version c(0.6, 0.7, 0.8)
OR2 = c(0.7,0.8)
p_init = c(1)

# scenarios = expand_grid(p0_e1,p0_e2,OR1,OR2,p_init)
# scenarios$scenario = 1:dim(scenarios)[1]
scenarios = expand_grid(p0_e1,p0_e2,OR1,OR2)
scenarios$scenario = 1:dim(scenarios)[1]
scenarios = expand_grid(scenarios,p_init)

# Probabilities treat group
scenarios$p1_e1 = (scenarios$OR1*scenarios$p0_e1/(1-scenarios$p0_e1))/(1+(scenarios$OR1*scenarios$p0_e1/(1-scenarios$p0_e1)))
scenarios$p1_e2 = (scenarios$OR2*scenarios$p0_e2/(1-scenarios$p0_e2))/(1+(scenarios$OR2*scenarios$p0_e2/(1-scenarios$p0_e2)))

# Calculate the correlation bounds
scenarios$min_corr0 = mapply(lower_corr,scenarios$p0_e1,scenarios$p0_e2)
scenarios$min_corr1 = mapply(lower_corr,scenarios$p1_e1,scenarios$p1_e2)
scenarios$max_corr0 = mapply(upper_corr,scenarios$p0_e1,scenarios$p0_e2)
scenarios$max_corr1 = mapply(upper_corr,scenarios$p1_e1,scenarios$p1_e2)

scenarios$min_corr = pmax(scenarios$min_corr0,scenarios$min_corr1)
scenarios$max_corr = pmin(scenarios$max_corr0,scenarios$max_corr1) 

corr = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8)
corr0 = c(0.2)

dataset = expand_grid(scenarios, corr, corr0) 
dataset = subset(dataset, dataset$corr < dataset$max_corr & dataset$corr > dataset$min_corr)

# save.image("C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/CBE_selection/results/scenarios.RData")

# Calculate CE Probabilities
dataset$p0_ce = mapply(prob_cbe, p_e1=dataset$p0_e1, p_e2=dataset$p0_e2, rho=dataset$corr0)
dataset$p1_ce = mapply(prob_cbe, p_e1=dataset$p1_e1, p_e2=dataset$p1_e2, rho=dataset$corr) 
dataset$OR_ce = (dataset$p1_ce/(1-dataset$p1_ce))/(dataset$p0_ce/(1-dataset$p0_ce)) 

# Calculate Sample size (total sample size, n=2*n0=2*n1)
alpha=0.05; beta=0.2
dataset$samplesize_e1 = mapply(samplesize_OR,p0=dataset$p0_e1, OR=dataset$OR1, alpha=alpha, beta=beta) 
dataset$samplesize_ce = mapply(samplesize_OR,p0=dataset$p0_ce, OR=dataset$OR_ce, alpha=alpha, beta=beta) 
dataset$samplesize_ce0 = mapply(samplesize_cbe,p0_e1=dataset$p0_e1,p0_e2=dataset$p0_e2,
                                eff_e1 = dataset$OR1,
                                effm_e1 = "or",
                                eff_e2 = dataset$OR2,
                                effm_e2 = "or",
                                effm_ce = "or",
                                rho=0, alpha=alpha, beta=beta)  
# dataset$ss_ratio = dataset$samplesize_e1/dataset$samplesize_ce
# dataset$ss_decision = ifelse(dataset$ss_ratio<1, "RE", "CE") 

# clean
rm(OR1,OR2,p0_e1,p0_e2,corr,corr0,p_init,scenarios)

# 

# Vector empirical powers and significance level
dataset$Test_Power_CE <- NA
dataset$Test_Power_RE <- NA

# blinded 
dataset$Test_Power_ES_SS <- NA
dataset$decision_ES_SS <- NA

#########################################
# General Settings Simulations
#########################################


library(doParallel)
# setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl) 

set.seed(156)


#########################################

# nsim: number of simulations
nsim = 100000
# nsim = 1000#test

# type i and ii errors
z.alpha <- qnorm(1-alpha,0,1)  
z.beta <-  qnorm(1-beta,0,1) 



#########################################
# Simulation under H1 
# - without sample size reassessment  
#########################################

set.seed(5532)

t0=Sys.time() 

#  COMPOSITE and RELEVANT ENDPOINT - Powers (for comparison)

for(i in 1:dim(dataset)[1]){
  dataset$Test_Power_CE[i] <- sum(replicate(nsim,f_OR(dataset$samplesize_e1[i]/2,dataset$p0_ce[i],dataset$p1_ce[i]))< - z.alpha)/nsim
  dataset$Test_Power_RE[i] <- sum(replicate(nsim,f_OR(dataset$samplesize_e1[i]/2,dataset$p0_e1[i],dataset$p1_e1[i]))< - z.alpha)/nsim
  print(i)
}

#  ENDPOINT SELECTION - Blinded approach 

for(i in 1:dim(dataset)[1]){
  aux <- rowSums(replicate(nsim,
                           f_dcorr(ss_arm=round(dataset$samplesize_e1[i]/2),
                                   p0_e1=dataset$p0_e1[i],OR1=dataset$OR1[i],
                                   p0_e2=dataset$p0_e2[i],OR2=dataset$OR2[i],
                                   p0_ce=dataset$p0_ce[i],p1_ce=dataset$p1_ce[i], 
                                   p_init=dataset$p_init[i],criteria="SS"))< c(-z.alpha,1))/nsim 
  
  
  dataset$Test_Power_ES_SS[i]<- aux[1]
  dataset$decision_ES_SS[i]<- 1-aux[2]
  print(i)
} 


t1=Sys.time()-t0   
# save.image("C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/CBE_selection/CODE_paper/simulations/results/results_diffcorr.RData") 

##################################################################################### 
#####################################################################################


# stop cluster
stopCluster(cl)


##################################################################################### 
#####################################################################################
load("C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/CBE_selection/CODE_paper/simulations/results/results_diffcorr.RData") 

# selecting only those scenarios without interim analysis
dataset$ss_decision = ifelse(dataset$decision_ES_SS<1, "RE", "CE")
dataset_H0False_ss = subset(dataset,dataset$p_init==1) 
# rm(dataset)

# 
summary(dataset_H0False_ss)
 
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

data_plot = subset(dataset_H0False_ss,dataset_H0False_ss$p_init==1.0) 
power_data_blinded_plot <- data.frame(Power=c(data_plot$Test_Power_ES_SS, 
                                              data_plot$Test_Power_CE,data_plot$Test_Power_RE),
                                 Design=c(rep("AD",length(data_plot$Test_Power_ES_SS)), 
                                          rep("CD",length(data_plot$Test_Power_CE)),
                                          rep("RD",length(data_plot$Test_Power_RE))),
                                 scenario=c(data_plot$scenario,
                                            data_plot$scenario, data_plot$scenario),
                                 corr=c(data_plot$corr,
                                        data_plot$corr, data_plot$corr)
)

for(i in 1:max(power_data_blinded_plot$scenario)){
  sub_leg=subset(dataset_H0False_ss,dataset_H0False_ss$scenario==i & dataset_H0False_ss$p_init==1.0) 
  sub=subset(power_data_blinded_plot,power_data_blinded_plot$scenario==i) 
  # univ=rbind(sub,sub)
  # univ$testind=c(rep("RE",dim(sub)[1]),rep("CE",dim(sub)[1]))
  # univ$test=ifelse(univ$testind=="RE",univ$Test_Power_RE,univ$Test_Power_CE)
  # univ$ss_decision = ifelse(dataset$decision_ES_SS<1, "RE", "CE") 
  # 
  # p[[enum]] <-ggplot(sub, aes(x=corr, y=Power))  +
  #   geom_point(size=2)+ ggtitle(paste("Scenario", sub_leg$scenario[it], "\n (p1,p2,OR1,OR2) \n=(", sub_leg$p0_e1[it],",",sub_leg$p0_e2[it],",",sub_leg$OR1[it],",",sub_leg$OR2[it],")"))+geom_point(size=2)  + labs(y = "Power (AD)", x="Correlation", color="Decision SS") + coord_cartesian(ylim = c(0.40, 1))+ geom_path()+ theme(plot.title = element_text(size=9),legend.position="bottom", legend.title = element_text(size = 6), legend.text = element_text(size = 6))
  # + theme(plot.title = element_text(size=9),legend.position = c(0.8, 0.2)) 
  
  #
  p[[enum]] <- ggplot(sub, aes(x=corr, y=Power, color=as.factor(Design)))+
    geom_point(size=2)+ ggtitle(paste("Scenario", sub_leg$scenario[it], "\n (p1,p2,OR1,OR2) \n=(", sub_leg$p0_e1[it],",",sub_leg$p0_e2[it],",",sub_leg$OR1[it],",",sub_leg$OR2[it],")"))+geom_point(size=2)  + labs(y = "Power", x="Correlation in the treatment arm", color="Design") + coord_cartesian(ylim = c(0.40, 1))+ geom_path()+ theme(plot.title = element_text(size=9),legend.position="bottom", legend.title = element_text(size = 6), legend.text = element_text(size = 6))+ geom_vline(xintercept = 0.2,linetype=2)
  
  enum=enum+1 
  # it <- it + dim(subset(dataset_H0False_ss,dataset_H0False_ss$scenario==i))
}
plot <- marrangeGrob(p,ncol=2,nrow=1,top=NULL)  
ggsave(file="results_diffcorr.pdf", plot, width = 210, height = 155, units = "mm") 


##################################################################################### 
##################################################################################### 

# Calculate the sample size with respect to $\rho^{(1)}$

# Study on how the sample size varies with respect to $\rho^{(1)}$. The black line represents the sample size computed by allowing different correlations, blue line represents the sample size assuming that the correlations are the same in both arms and equal to  $\rho^{(0)}$, and  the red line is the sample size assuming that the correlations are the same in both arms and equal to  $(\rho^{(0)}+\rho^{(1)})/2$.

p0_e1=0.1;p0_e2=0.25;p1_e1=0.06;p1_e2=0.19;alpha=0.05;beta=0.2;rho0=0.2;step=5
rho1=seq(0,max,length.out=step)
max=min(upper_corr(p_e1=p0_e1, p_e2=p0_e2),upper_corr(p_e1=p1_e1, p_e2=p1_e2))

sss=Vectorize(
  function(rho1,rho0=0.2,p0_e1=0.1,p0_e2=0.25,p1_e1=0.06,p1_e2=0.21,alpha=0.05,beta=0.2){
    
    # max=min(upper_corr(p_e1=p0_e1, p_e2=p0_e2),upper_corr(p_e1=p1_e1, p_e2=p1_e2))
    # rho1=seq(0,max,length.out=step)
    
    p0_ce=prob_cbe(p_e1=p0_e1, p_e2=p0_e2, rho=rho0)
    p1_ce=prob_cbe(p_e1=p1_e1, p_e2=p1_e2, rho=rho1)
    OR_ce=(p1_ce/(1-p1_ce))/(p0_ce/(1-p0_ce)) 
    ss=samplesize_OR(p0=p0_ce, OR=OR_ce, alpha=alpha, beta=beta)
    
    p0_ce_av=prob_cbe(p_e1=p0_e1, p_e2=p0_e2, rho=(rho0+rho1)/2)
    p1_ce_av=prob_cbe(p_e1=p1_e1, p_e2=p1_e2, rho=(rho0+rho1)/2)
    OR_ce_av=(p1_ce_av/(1-p1_ce_av))/(p0_ce_av/(1-p0_ce_av)) 
    ss_av=samplesize_OR(p0=p0_ce_av, OR=OR_ce_av, alpha=alpha, beta=beta)
    
    p0_ce_eq=prob_cbe(p_e1=p0_e1, p_e2=p0_e2, rho=rho0)
    p1_ce_eq=prob_cbe(p_e1=p1_e1, p_e2=p1_e2, rho=rho0)
    OR_ce_eq=(p1_ce_eq/(1-p1_ce_eq))/(p0_ce_eq/(1-p0_ce_eq)) 
    ss_eq=samplesize_OR(p0=p0_ce_eq, OR=OR_ce_eq, alpha=alpha, beta=beta)
    
    # c(rho0,rho1,ss,ss_av,ss_eq)
    
    list(rho0,rho1,ss,ss_av,ss_eq,p0_ce,p1_ce,OR_ce)
    
  },vectorize.args = "rho1")

sss(rho1=0.5)
(ss_example = sss(rho1=seq(0.1,0.3,length.out=10)))


# Plots
plot(ss_example[2,],ss_example[3,],type="l",ylab="Sample sizes",xlab="Correlation in arm 1")#,xlim=c(0,0.5) 
lines(ss_example[2,],ss_example[4,],col="red")
lines(ss_example[2,],ss_example[5,],col="blue") 
