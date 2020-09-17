
##################################################################################
# Research project Vienna - Bcn
# Endpoint selection and sample size reassessment for composite binary endpoints
# Simulation study 
##################################################################################

setwd("C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/CBE_selection")
source('C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/CBE_selection/sim_functions.R')
source('C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/CBE_selection/sim_functions_ce.R')

# setwd("C:/Users/Marta/Dropbox/C5/Scripts/GitKraken/CBE_selection")
# source('C:/Users/Marta/Dropbox/C5/Scripts/GitKraken/CBE_selection/sim_functions.R')
# source('C:/Users/Marta/Dropbox/C5/Scripts/GitKraken/CBE_selection/sim_functions_ce.R')

#########################################
# Preamble
#########################################

library(tidyverse)
library(tidyr)
library(plyr)
library(devtools) 
install_github("CompARE-Composite/CompARE-package")
library(CompARE)

#########################################
# Define the set of scenarios
#########################################

# Scenarios
p0_e1 = c(0.1, 0.2)
p0_e2 = c(0.1, 0.2, 0.3)
OR1 = c(0.6, 0.7, 0.8, 0.9)
OR2 = c(0.65, 0.7, 0.8)

scenarios = expand_grid(p0_e1,p0_e2,OR1,OR2)
scenarios$scenario = 1:dim(scenarios)[1]

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

corr = c(0, 0.2, 0.4, 0.6, 0.8)

dataset = expand_grid(scenarios, corr) 
dataset = subset(dataset, dataset$corr < dataset$max_corr & dataset$corr > dataset$min_corr)

# Calculate CE Probabilities
dataset$p0_ce = mapply(prob_cbe, p_e1=dataset$p0_e1, p_e2=dataset$p0_e2, rho=dataset$corr)
dataset$p1_ce = mapply(prob_cbe, p_e1=dataset$p1_e1, p_e2=dataset$p1_e2, rho=dataset$corr) 
dataset$OR_ce = (dataset$p1_ce/(1-dataset$p1_ce))/(dataset$p0_ce/(1-dataset$p0_ce))

# Calculate ARE
dataset$ARE = mapply(ARE_cbe, p0_e1=dataset$p0_e1, p0_e2=dataset$p0_e2, eff_e1=dataset$OR1, eff_e2=dataset$OR2, rho=dataset$corr)
dataset$decision = ifelse(dataset$ARE<1, "RE", "CE")

rm(OR1,OR2,p0_e1,p0_e2,corr,scenarios)

#########################################
# Simulations
#########################################

set.seed(4123)

# nsim: number of simulations
nsim = 100000

# sample size per arm n0=n1
n0 = 500

# type i and ii errors
alpha=0.025; beta=0.2
z.alpha <- qnorm(1-alpha,0,1)  
z.beta <-  qnorm(1-beta,0,1) 


t0=Sys.time()

dataset$Test_Power_CE = 0
dataset$Test_Power_RE = 0
dataset$Test_Power_ES = 0 

dataset$Test_Reject_CE = 0
dataset$Test_Reject_RE = 0
dataset$Test_Reject_ES = 0 

#########################################

for(i in 1:dim(dataset)[1]){
  dataset$Test_Power_CE[i] <- sum(replicate(nsim,f_OR(n0,dataset$p0_ce[i],dataset$p1_ce[i]))< - z.alpha)/nsim
  dataset$Test_Power_RE[i] <- sum(replicate(nsim,f_OR(n0,dataset$p0_e1[i],dataset$p1_e1[i]))< - z.alpha)/nsim
  print(i)
}

save.image("C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/CBE_selection/results/results.RData")

######################################### 

for(i in 1:dim(dataset)[1]){ 
  dataset$Test_Power_ES[i] <- tryCatch(sum(replicate(nsim,
                                            f_ES(samplesize=2*n0,
                                                 p0_e1=dataset$p0_e1[i],p1_e1=dataset$p1_e1[i],
                                                 OR1=dataset$OR1[i],
                                                 p0_e2=dataset$p0_e2[i],p1_e2=dataset$p1_e2[i],
                                                 OR2=dataset$OR2[i],
                                                 p0_ce=dataset$p0_ce[i],p1_ce=dataset$p1_ce[i],
                                                 OR_ce=dataset$OR_ce[i],
                                                 upp=dataset$max_corr[i],low=dataset$min_corr[i]))< - z.alpha)/nsim,
                                       error=function(e){NA})
  print(i)
}

t1=Sys.time()-t0   
save.image("C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/CBE_selection/results/results_decision.RData")
# save.image("C:/Users/Marta/Dropbox/C5/Scripts/GitKraken/CBE_selection/results_decision.RData")


#########################################

for(i in 1:dim(dataset)[1]){
  dataset$Test_Reject_CE[i] <- sum(replicate(nsim,f_OR(n0,dataset$p0_ce[i],dataset$p0_ce[i]))< - z.alpha)/nsim
  dataset$Test_Reject_RE[i] <- sum(replicate(nsim,f_OR(n0,dataset$p0_e1[i],dataset$p0_e1[i]))< - z.alpha)/nsim
  print(i)
}

save.image("C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/CBE_selection/results/results_H0.RData")

######################################### 

for(i in 1:dim(dataset)[1]){ 
  dataset$Test_Reject_ES[i] <- tryCatch(sum(replicate(nsim,
                                             f_ES(samplesize=2*n0,
                                                  p0_e1=dataset$p0_e1[i],p1_e1=dataset$p0_e1[i],
                                                  OR1=dataset$OR1[i],
                                                  p0_e2=dataset$p0_e2[i],p1_e2=dataset$p0_e2[i],
                                                  OR2=dataset$OR2[i],
                                                  p0_ce=dataset$p0_ce[i],p1_ce=dataset$p0_ce[i],
                                                  OR_ce=dataset$OR_ce[i],
                                                  upp=dataset$max_corr[i],low=dataset$min_corr[i]))< - z.alpha)/nsim,
                                        error=function(e){NA})
  print(i)
}

t1=Sys.time()-t0   
save.image("C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/CBE_selection/results/results_decision_H0.RData")
# save.image("C:/Users/Marta/Dropbox/C5/Scripts/GitKraken/CBE_selection/results_decision_H0.RData")

######################################### 
#########################################  

