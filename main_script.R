
##################################################################################
# Research project Vienna - Bcn
# Endpoint selection and sample size reassessment for composite binary endpoints
# Simulation study 
##################################################################################

library(tidyverse)
library(tidyr)
library(plyr)
library(devtools) 
install_github("CompARE-Composite/CompARE-package")
library(CompARE)

# Scenarios
p0_e1 = c(0.1, 0.3)
p0_e2 = c(0.2, 0.4)
OR1 = c(0.7, 0.8, 0.9)
OR2 = c(0.6, 0.7, 0.8)

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



##################################################################################

# check OR transformation
OR.function<- function(p0, p1){
  OR<- (p1/(1-p1))/(p0/(1-p0))
  return(OR)
}

scenarios$OR = mapply(OR.function,scenarios$p0_e1,scenarios$p1_e1) 
(scenarios$OR == scenarios$OR1)
(scenarios$OR - scenarios$OR1)
