
##################################################################################
# Research project Vienna - Bcn
# Endpoint selection and sample size reassessment for composite binary endpoints
# Simulation study 
##################################################################################

rm(list = ls())

# setwd("H:/Code_PAPER_biost")
# source('H:/Code_PAPER_biost/aux_functions.R') 
# source('H:/Code_PAPER_biost/eselect/R/eselect.R')
# source('H:/Code_PAPER_biost/eselect/R/eselectsim.R') 
# source('H:/Code_PAPER_biost/eselect/R/eselect_ub.R')
# source('H:/Code_PAPER_biost/eselect/R/eselectsim_ub.R') 

setwd("~/Code_PAPER_biost")
source('~/Code_PAPER_biost/aux_functions.R') 
source('~/Code_PAPER_biost/eselect/R/eselect.R')
source('~/Code_PAPER_biost/eselect/R/eselectsim.R') 
source('~/Code_PAPER_biost/eselect/R/eselect_ub.R')
source('~/Code_PAPER_biost/eselect/R/eselectsim_ub.R') 


# 
# setwd("C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/CBE_selection")
# source('C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/CBE_selection/aux_functions.R') 
# source('C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/CBE_selection/eselect/R/eselect.R')
# source('C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/CBE_selection/eselect/R/eselectsim.R') 
# source('C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/CBE_selection/eselect/R/eselect_ub.R')
# source('C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/CBE_selection/eselect/R/eselectsim_ub.R') 

##################################################################################
# Differences with respect to previous versions
##################################################################################

# - consider sample size reassessment with initial sample size cbe
# - consider scenarios without sample size reassessment (and initial sample size e1)

#########################################
# Preamble
#########################################

library(doParallel)
# setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl) 

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

#########################################
# Define the set of scenarios
#########################################

# Scenarios
p0_e1 = c(0.1, 0.2)
p0_e2 = c(0.1, 0.25)
OR1 = c(0.6, 0.8) #old_version c(0.6, 0.7, 0.8)
OR2 = c(0.75, 0.8)
p_init = c(0.5, 1)

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

# for mild cases
scenarios$min_corr0_1 = mapply(lower_corr,scenarios$p0_e1,scenarios$p1_e2)
scenarios$min_corr1_0 = mapply(lower_corr,scenarios$p1_e1,scenarios$p0_e2)
scenarios$max_corr0_1 = mapply(upper_corr,scenarios$p0_e1,scenarios$p1_e2)
scenarios$max_corr1_0 = mapply(upper_corr,scenarios$p1_e1,scenarios$p0_e2)
# 
scenarios$min_corr = pmax(scenarios$min_corr0,scenarios$min_corr1,scenarios$min_corr1_0,scenarios$min_corr0_1)
scenarios$max_corr = pmin(scenarios$max_corr0,scenarios$max_corr1,scenarios$max_corr1_0,scenarios$max_corr0_1) 

corr = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)

dataset = expand_grid(scenarios, corr) 
dataset = subset(dataset, dataset$corr < dataset$max_corr & dataset$corr > dataset$min_corr)

# save.image("C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/CBE_selection/results/scenarios.RData")

# Calculate CE Probabilities
dataset$p0_ce = mapply(prob_cbe, p_e1=dataset$p0_e1, p_e2=dataset$p0_e2, rho=dataset$corr)
dataset$p1_ce = mapply(prob_cbe, p_e1=dataset$p1_e1, p_e2=dataset$p1_e2, rho=dataset$corr) 
dataset$OR_ce = (dataset$p1_ce/(1-dataset$p1_ce))/(dataset$p0_ce/(1-dataset$p0_ce))

# Calculate ARE
dataset$ARE = mapply(ARE_cbe, p0_e1=dataset$p0_e1, p0_e2=dataset$p0_e2, eff_e1=dataset$OR1, eff_e2=dataset$OR2, rho=dataset$corr)
dataset$decision = ifelse(dataset$ARE<1, "RE", "CE") 

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

dataset$ss_ratio = dataset$samplesize_e1/dataset$samplesize_ce
dataset$ss_decision = ifelse(dataset$ss_ratio<1, "RE", "CE") 

# clean
rm(OR1,OR2,p0_e1,p0_e2,corr,p_init,scenarios)

# Vector empirical powers and significance level
dataset$Test_Power_CE <- NA
dataset$Test_Power_RE <- NA

# blinded 
dataset$Test_Power_ES_SS <- NA
dataset$decision_ES_SS <- NA

# unblinded 
dataset$Test_Power_ES_ubSS <- NA 
dataset$decision_ES_ubSS  <- NA 

# Empirical alpha
dataset$Test_Reject_CE <- NA
dataset$Test_Reject_RE <- NA
dataset$Test_Reject_ES <- NA
dataset$Test_Reject_ES_ub <- NA

#########################################
# General Settings Simulations
#########################################

set.seed(4325)

# nsim: number of simulations
nsim = 100000
# nsim = 1000#test

# type i and ii errors
z.alpha <- qnorm(1-alpha,0,1)  
z.beta <-  qnorm(1-beta,0,1) 


#########################################
# Simulation under H1
#########################################

t0=Sys.time()   

#########################################
# Simulation under H1 
# Simulating that one endpoint (endpoint2) has no effect 
######################################### 


set.seed(23423)

t0=Sys.time() 

# prob_cbe(p_e1=dataset$p0_e1[i], p_e2=dataset$p0_e2[i], rho=dataset$corr[i])
# prob_cbe(p_e1=dataset$p1_e1[i], p_e2=dataset$p0_e2[i], rho=dataset$corr[i])
# dataset$p0_ce = mapply(prob_cbe, p_e1=dataset$p0_e1, p_e2=dataset$p0_e2, rho=dataset$corr)
# dataset$p1_ce = mapply(prob_cbe, p_e1=dataset$p1_e1, p_e2=dataset$p1_e2, rho=dataset$corr) 

#  COMPOSITE and RELEVANT ENDPOINT - Powers (for comparison)

for(i in 1:dim(dataset)[1]){
  dataset$Test_Power_CE[i] <- sum(replicate(nsim,f_OR(dataset$samplesize_e1[i]/2,
                                                      prob_cbe(p_e1=dataset$p0_e1[i], p_e2=dataset$p0_e2[i], 
                                                               rho=dataset$corr[i]),
                                                      prob_cbe(p_e1=dataset$p1_e1[i], p_e2=dataset$p0_e2[i], 
                                                               rho=dataset$corr[i])))< - z.alpha)/nsim
  dataset$Test_Power_RE[i] <- sum(replicate(nsim,f_OR(dataset$samplesize_e1[i]/2,dataset$p0_e1[i],dataset$p1_e1[i]))< - z.alpha)/nsim
  print(i)
}

#  ENDPOINT SELECTION - Blinded approach 

for(i in 1:dim(dataset)[1]){
  aux <- rowSums(replicate(nsim,eselectsim(ss_arm=round(dataset$samplesize_e1[i]/2),
                                           p0_e1=dataset$p0_e1[i],OR1=dataset$OR1[i],
                                           p0_e2=dataset$p0_e2[i],OR2=dataset$OR2[i],
                                           p0_ce=dataset$p0_ce[i],p_init=dataset$p_init[i],
                                           criteria="SS",H0_e1=FALSE,H0_e2=TRUE))< c(-z.alpha,1))/nsim 
  
  
  dataset$Test_Power_ES_SS[i]<- aux[1]
  dataset$decision_ES_SS[i]<- 1-aux[2]
  print(i)
}

#  ENDPOINT SELECTION - Unblinded approach  

for(i in 1:dim(dataset)[1]){
  aux <- rowSums(replicate(nsim,eselectsim_ub(ss_arm=round(dataset$samplesize_e1[i]/2),
                                              p0_e1=dataset$p0_e1[i],OR1=dataset$OR1[i],
                                              p0_e2=dataset$p0_e2[i],OR2=dataset$OR2[i],
                                              p0_ce=dataset$p0_ce[i],p_init=dataset$p_init[i],
                                              criteria="SS",H0_e1=FALSE,H0_e2=TRUE))< c(-z.alpha,1))/nsim 
  
  
  dataset$Test_Power_ES_ubSS[i]<- aux[1]
  dataset$decision_ES_ubSS[i]<- 1-aux[2]
  print(i)
}

t1=Sys.time()-t0  
# save.image("H:/Code_PAPER_biost/results/results_H0True_e2.RData") 
save.image("~/Code_PAPER_biost/results/results_H0True_e2.RData") 


#########################################
# Simulation under H1 
# Simulating that one endpoint (endpoint1) has no effect 
######################################### 

t0=Sys.time()  

#  COMPOSITE and RELEVANT ENDPOINT - Powers (for comparison)

for(i in 1:dim(dataset)[1]){
  dataset$Test_Power_CE[i] <- sum(replicate(nsim,f_OR(dataset$samplesize_e1[i]/2,
                                                      prob_cbe(p_e1=dataset$p0_e1[i], p_e2=dataset$p0_e2[i], 
                                                               rho=dataset$corr[i]),
                                                      prob_cbe(p_e1=dataset$p0_e1[i], p_e2=dataset$p1_e2[i], 
                                                               rho=dataset$corr[i])))< - z.alpha)/nsim
  dataset$Test_Power_RE[i] <- sum(replicate(nsim,f_OR(dataset$samplesize_e1[i]/2,dataset$p0_e1[i],dataset$p0_e1[i]))< - z.alpha)/nsim
  print(i)
}

#  ENDPOINT SELECTION - Blinded approach 

for(i in 1:dim(dataset)[1]){
  aux <- rowSums(replicate(nsim,eselectsim(ss_arm=round(dataset$samplesize_e1[i]/2),
                                           p0_e1=dataset$p0_e1[i],OR1=dataset$OR1[i],
                                           p0_e2=dataset$p0_e2[i],OR2=dataset$OR2[i],
                                           p0_ce=dataset$p0_ce[i],p_init=dataset$p_init[i],
                                           criteria="SS",H0_e1=TRUE,H0_e2=FALSE))< c(-z.alpha,1))/nsim 
  
  
  dataset$Test_Power_ES_SS[i]<- aux[1]
  dataset$decision_ES_SS[i]<- 1-aux[2]
  print(i)
}

#  ENDPOINT SELECTION - Unblinded approach  

for(i in 1:dim(dataset)[1]){
  aux <- rowSums(replicate(nsim,eselectsim_ub(ss_arm=round(dataset$samplesize_e1[i]/2),
                                              p0_e1=dataset$p0_e1[i],OR1=dataset$OR1[i],
                                              p0_e2=dataset$p0_e2[i],OR2=dataset$OR2[i],
                                              p0_ce=dataset$p0_ce[i],p_init=dataset$p_init[i],
                                              criteria="SS",H0_e1=TRUE,H0_e2=FALSE))< c(-z.alpha,1))/nsim 
  
  
  dataset$Test_Power_ES_ubSS[i]<- aux[1]
  dataset$decision_ES_ubSS[i]<- 1-aux[2]
  print(i)
}

t1=Sys.time()-t0  
# save.image("H:/Code_PAPER_biost/results/results_H0True_e1.RData") 
save.image("~/Code_PAPER_biost/results/results_H0True_e1.RData") 


##################################################################################### 
#####################################################################################


# stop cluster
stopCluster(cl)
