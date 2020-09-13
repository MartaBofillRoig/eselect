
##################################################################################
# Research project Vienna - Bcn
# Endpoint selection and sample size reassessment for composite binary endpoints
# Simulation study 
##################################################################################

setwd("C:/Users/Marta/Dropbox/C5/Scripts/GitKraken/CBE_selection")
source('C:/Users/Marta/Dropbox/C5/Scripts/GitKraken/CBE_selection/sim_functions.R')

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

scenarios1 = scenarios; scenarios2 = scenarios; scenarios3 = scenarios

scenarios3$corr = (scenarios3$max_corr)/3
scenarios2$corr = (scenarios2$max_corr)/2
scenarios1$corr = scenarios1$max_corr
scenarios$corr = 0

dataset = rbind(scenarios,scenarios1,scenarios2,scenarios3)

# Calculate CE Probabilities
dataset$p0_ce = mapply(prob_cbe, p_e1=dataset$p0_e1, p_e2=dataset$p0_e2, rho=dataset$corr)
dataset$p1_ce = mapply(prob_cbe, p_e1=dataset$p1_e1, p_e2=dataset$p1_e2, rho=dataset$corr) 

# Calculate ARE
dataset$ARE = mapply(ARE_cbe, p0_e1=dataset$p0_e1, p0_e2=dataset$p0_e2, eff_e1=dataset$OR1, eff_e2=dataset$OR2, rho=dataset$corr)

dataset$decision = ifelse(dataset$ARE<1, "RE", "CE")

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

#########################################
t0=Sys.time()

dataset$Test_Reject_CE = 0
dataset$Test_Reject_RE = 0

for(i in 1:dim(dataset)[1]){
  dataset$Test_Reject_CE[i] <- sum(replicate(nsim,f_OR(n0,dataset$p0_ce[i],dataset$p1_ce[i]))< - z.alpha)/nsim
  dataset$Test_Reject_RE[i] <- sum(replicate(nsim,f_OR(n0,dataset$p0_e1[i],dataset$p1_e1[i]))< - z.alpha)/nsim
  print(i)
}

t1=Sys.time()-t0 
rm(alpha,beta,i,nsim,z.alpha,z.beta,OR1,OR2,p0_ce,p0_e1,p0_e2,f_OR,scenarios,scenarios1,scenarios2,scenarios3)
save.image("C:/Users/Marta/Dropbox/C5/Scripts/GitKraken/CBE_selection/results.RData")

#########################################
# Results
#########################################

dataset$diff_powers = dataset$Test_Reject_CE - dataset$Test_Reject_RE
dataset$gain_power = ifelse(dataset$diff_powers>0, "CE", "RE")

summary(dataset)

(dataset$gain_power == dataset$decision)
