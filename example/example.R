
##################################################################################
# Research project Vienna - Bcn
# Endpoint selection and sample size reassessment for composite binary endpoints
# Simulation study 
##################################################################################

rm(list = ls())

library(CompAREdesign)

setwd("C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/CBE_selection")
source('C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/CBE_selection/aux_functions.R') 
source('C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/CBE_selection/eselect/R/eselect.R')
source('C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/CBE_selection/eselect/R/eselectsim.R') 
source('C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/CBE_selection/eselect/R/eselect_ub.R')
source('C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/CBE_selection/eselect/R/eselectsim_ub.R') 


######################################
# Probabilities in the control group for E1, E2 and E3
p0_e1 = 0.15
p0_e2 = 0.50
p0_e3 = 0.23


######################################
# Overlap between components E2 and E3 
# Bounds for the relative overlap
ro0_b = c(max(0,p0_e2+p0_e3-1),min(p0_e2/p0_e3,p0_e3/p0_e2))

ro0_e23 = 0.4
p0_e23 = (ro0_e23*(p0_e2+p0_e3))/(1+ro0_e23)

corr0 = (p0_e23 - p0_e2*p0_e3)/sqrt(p0_e2*(1-p0_e2)*p0_e3*(1-p0_e3))
corr0

lower_corr(p_e1=p0_e2, p_e2=p0_e3)
upper_corr(p_e1=p0_e2, p_e2=p0_e3)

prob_cbe(p_e1=p0_e2, p_e2=p0_e3, rho=corr0)

######################################
#  Assuming correlation =0

p0_ce1 = prob_cbe(p_e1=p0_e2, p_e2=p0_e3, rho=0)

# Correlation bounds adding third component:
lower_corr(p_e1=p0_ce1, p_e2=p0_e1)
upper_corr(p_e1=p0_ce1, p_e2=p0_e1)




# 
prob_cbe(p_e1=p0_e1, p_e2=prob_cbe(p_e1=p0_e2, p_e2=p0_e3, rho=0), rho=0)

samplesize_cbe(
  p0_e1=p0_e1,
  p0_e2=prob_cbe(p_e1=p0_e2, p_e2=p0_e3, rho=0),
  eff_e1=0.66,
  effm_e1="or",
  eff_e2=0.47,
  effm_e2="or",
  effm_ce = "or",
  rho=c(0,0.1,0.2,0.3,0.4,0.5),
  alpha = 0.05,
  beta = 0.2,
  unpooled = TRUE
)

plot_ce(prob_cbe(p_e1=p0_e1, p_e2=prob_cbe(p_e1=p0_e2, p_e2=p0_e3, rho=0), rho=c(0,0.1,0.2,0.3)))



# Comparison e1 vs ce

# db <- table 2x2 results
# eselect(db,p0_e1,OR1=0.7,p0_e2,OR2=0.8,criteria="SS",alpha=0.05,beta=0.2)

eselectsim(ss_arm=500,p0_e1=p0_e1,OR1=0.66,p0_e2=p0_e2,OR2=0.47,
           p0_ce=prob_cbe(p_e1=p0_e2, p_e2=p0_e3, rho=0),
           p_init=1,criteria="SS",SS_r=FALSE,alpha=0.05,beta=0.2)




