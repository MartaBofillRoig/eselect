##################################################################################
# Research project Vienna - Bcn
# Endpoint selection and sample size reassessment for composite binary endpoints
# Example supp material
##################################################################################

rm(list = ls())

library(CompAREdesign)
library(tidyverse)
library(gridExtra)

setwd("C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/CBE_selection")
source('C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/CBE_selection/CODE_paper/simulations/aux_functions.R')
source('C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/CBE_selection/eselect/R/eselect.R')
source('C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/CBE_selection/eselect/R/eselectsim.R') 
source('C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/CBE_selection/eselect/R/eselect_ub.R')
source('C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/CBE_selection/eselect/R/eselectsim_ub.R') 


######################################
# EXAMPLE DATA
######################################

# pooled estimates
# Values in TAXUS V:
p0_e1 = 0.173
p0_e2 = 0.055

p1_e1 = 0.121; 
p1_e2 = 0.057;

n1 = 569
n0 = 576
n = n0+n1

p1 = (p0_e1*n0 + p1_e1*n1)/n
p2 = (p0_e2*n0 + p1_e2*n1)/n
p_ce = (0.203*n0 + 0.146*n1)/n

OR1 = 0.7
OR2 = 0.9

x11 = round((p1+p2-p_ce)*n)
x12 = round((p1)*n-x11)
x21 = round((p2)*n- x11)
x22 = round((1-p_ce)*n)

data = matrix(c(x11,x12,x21,x22), nrow = 2 , ncol = 2, byrow = F)

# eselect(db=data,p0_e1=0.18,OR1=0.65,p0_e2=0.05,OR2=0.9,criteria="SS",alpha=0.05,beta=0.2) 
eselect(db=data,p0_e1=0.18,OR1=0.70,p0_e2=0.05,OR2=0.9,criteria="SS",alpha=0.05,beta=0.2) 


# Simulations considering different correlation values
# eselect
set.seed(5234)
rho=c(0,0.1,0.2,0.3,0.4,0.5) 
(p0_ce = prob_cbe(p_e1=0.18, p_e2=0.05, rho=rho))

eselectsim(ss_arm=n/2, p0_e1=0.18, OR1=OR1, p0_e2=0.05, OR2=0.90, p0_ce=p0_ce[1], p_init = 1, alpha = 0.05, beta = 0.2)

eselectsim(ss_arm=n/2, p0_e1=0.18, OR1=OR1, p0_e2=0.05, OR2=0.90, p0_ce=p0_ce[3], p_init = 1, alpha = 0.05, beta = 0.2)

eselectsim(ss_arm=n/2, p0_e1=0.18, OR1=OR1, p0_e2=0.05, OR2=0.90, p0_ce=p0_ce[4], p_init = 1, alpha = 0.05, beta = 0.2)

eselectsim(ss_arm=n/2, p0_e1=0.18, OR1=OR1, p0_e2=0.05, OR2=0.80, p0_ce=p0_ce[1], p_init = 1, alpha = 0.05, beta = 0.2)

eselectsim(ss_arm=n/2, p0_e1=0.18, OR1=OR1, p0_e2=0.05, OR2=0.80, p0_ce=p0_ce[3], p_init = 1, alpha = 0.05, beta = 0.2)

eselectsim(ss_arm=n/2, p0_e1=0.18, OR1=OR1, p0_e2=0.05, OR2=0.80, p0_ce=p0_ce[4], p_init = 1, alpha = 0.05, beta = 0.2)



alpha =0.05
beta =0.2

samplesize_OR(p0=p0_e1, OR=OR1, alpha=alpha, beta=beta)

mapply(samplesize_cbe, p0_e1=p0_e1,
       p0_e2=p0_e2,
       eff_e1=OR1,
       effm_e1="or",
       eff_e2=0.80,
       effm_e2="or",
       effm_ce = "or",
       rho=c(0,0.1,0.2,0.3,0.4,0.5),alpha=alpha, beta=beta) 


######################################
# SAMPLE SIZE CALCULATIONS
######################################
# Probabilities in the control group for E1, E2 and E3

# Values in TAXUS V:
p0_e1 = 0.173
p0_e2 = 0.055

p1_e1 = 0.121; 
p1_e2 = 0.057;

a = c(0.057, 0.050, 0.045, 0.040)
OR1 = OR_function(p0=p0_e1,p1=p1_e1)
OR_function(p0=p0_e2,p1=a)
OR2 = 0.75

alpha =0.05
beta =0.2

n1 = samplesize_OR(p0=p0_e1, OR=OR1, alpha=alpha, beta=beta)

# Correlation bounds

lower_corr(p_e1=p0_e1, p_e2=p0_e2)
upper_corr(p_e1=p0_e1, p_e2=p0_e2)

# Composite endpoint

p0_ce = prob_cbe(p_e1=p0_e1, p_e2=p0_e2, rho=0)
OR_mace = effectsize_cbe(p0_e1=p0_e1, p0_e2=p0_e2,  
                         eff_e1=OR1,
                         effm_e1 = "or",
                         eff_e2=OR2,
                         effm_e2 = "or",
                         effm_ce = "or",
                         rho=0
)$`Effect CE`

# Sample size

rho=c(0,0.1,0.2,0.3,0.4,0.5) 

ss_cbe <- mapply(samplesize_cbe, p0_e1=p0_e1,
                 p0_e2=p0_e2,
                 eff_e1=OR1,
                 effm_e1="or",
                 eff_e2=0.75,
                 effm_e2="or",
                 effm_ce = "or",
                 rho=c(0,0.1,0.2,0.3,0.4,0.5),alpha=alpha, beta=beta) 

list(rho=rho,samplesize=ss_cbe)

# 
data_plot_simpl = data.frame(ss=c(ss_cbe, 
                                  rep(n1,length(ss_cbe)) 
),
rho=rep(rho,2),
Design=c(rep("CD",length(ss_cbe)),rep("RD",length(ss_cbe))
) 
)

data_plot_simpl$Design = factor(data_plot_simpl$Design,levels=c("CD","RD")) 

ggplot(data_plot_simpl,aes(x = rho, y = ss, color=Design#, linetype = Correlation
)) + geom_point() +geom_line(size=0.8) + xlab("Correlation") + ylab("Sample size")  + scale_color_manual(breaks = c("CD","RD"), values=c("#D9717D", "#4DB6D0")) 
ggsave("example_s_ex2.pdf",width = 100, height = 100, units = "mm")


# eselect
set.seed(5234)
(p0_ce = prob_cbe(p_e1=p0_e1, p_e2=p0_e2, rho=rho))

eselectsim( ss_arm=n1/2, p0_e1=p0_e1, OR1=OR1, p0_e2=p0_e2, OR2=OR2, p0_ce=p0_ce[1], p_init = 1, alpha = 0.05, beta = 0.2)

eselectsim(ss_arm=n1/2, p0_e1=p0_e1, OR1=OR1, p0_e2=p0_e2, OR2=OR2, p0_ce=p0_ce[3], p_init = 1, alpha = 0.05, beta = 0.2)

eselectsim(ss_arm=n1/2, p0_e1=p0_e1, OR1=OR1, p0_e2=p0_e2, OR2=OR2, p0_ce=p0_ce[4], p_init = 1, alpha = 0.05, beta = 0.2)

eselectsim( ss_arm=n1/2, p0_e1=p0_e1, OR1=OR1, p0_e2=p0_e2, OR2=0.80, p0_ce=p0_ce[1], p_init = 1, alpha = 0.05, beta = 0.2)

eselectsim(ss_arm=n1/2, p0_e1=p0_e1, OR1=OR1, p0_e2=p0_e2, OR2=0.80, p0_ce=p0_ce[3], p_init = 1, alpha = 0.05, beta = 0.2)

eselectsim(ss_arm=n1/2, p0_e1=p0_e1, OR1=OR1, p0_e2=p0_e2, OR2=0.80, p0_ce=p0_ce[4], p_init = 1, alpha = 0.05, beta = 0.2)


