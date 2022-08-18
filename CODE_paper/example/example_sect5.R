##################################################################################
# Code for the paper
# ``Adaptive clinical trial designs with blinded selection of binary composite endpoints and sample size reassessment''
# Motivating example - Sect 5
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
# Values based on Table 1 in Boehm et al.
# (Boehm, M., Niewczas, J., Herkner, H., Koenig, F., Kratochwill, K., Rutherford, P., Aufricht, C., & Vychytil, A. (2019). Composite Outcome Improves Feasibility of Clinical Trials in Peritoneal Dialysis. Peritoneal Dialysis International: Journal of the International Society for Peritoneal Dialysis, 39(5), 479–485. https://doi.org/10.3747/pdi.2018.00214)

# Probabilities in the control group for E1, E2 and E3
p0_e1 = 0.15 #E1: technical failure
p0_e2 = 0.50 #E2: peritonitis
p0_e3 = 0.23 #E3: peritoneal membrane deterioration

# Odds ratios
OR1=0.66
OR2=0.54
OR3=0.64

alpha=0.05
beta=0.2

# For the example in the paper, we consider the composite endpoint of peritonitis and peritoneal membrane deterioration (E2 U E3) and technical failure (E1).
# Note that we grouped peritonitis and peritoneal events together for the sake of illustration

# As in Boehm et al, peritonitis and peritoneal membrane deterioration is obtained assuming correlation =0
p0_ce1 = prob_cbe(p_e1=p0_e2, p_e2=p0_e3, rho=0)
p0_ce1

# Correlation bounds adding third component:
lower_corr(p_e1=p0_ce1, p_e2=p0_e1)
upper_corr(p_e1=p0_ce1, p_e2=p0_e1)

# We consider the relevant endpoint to be the endpoint of peritonitis and peritoneal membrane deterioration  (E2 U E3)

######################################
# Relevant endpoint
p0_ce1 = prob_cbe(p_e1=p0_e2, p_e2=p0_e3, rho=0) 
# Compute the odds ratio
OR_ce1 = effectsize_cbe(p0_e1=p0_e2, p0_e2=p0_e3,  
                        eff_e1=OR2,
                        effm_e1 = "or",
                        eff_e2=OR3,
                        effm_e2 = "or",
                        effm_ce = "or",
                        rho=0
)$`Effect CE`

# Compute sample size
n1 = samplesize_OR(p0=p0_ce1, OR=OR_ce1, alpha=alpha, beta=beta)

######################################
# Composite endpoint

# Calculate the correlation bounds between the components
lower_corr(p_e1=p0_ce1, p_e2=p0_e3)
upper_corr(p_e1=p0_ce1, p_e2=p0_e3)


prob_cbe(p_e1=p0_ce1, p_e2=p0_e3, rho=0) 

OR_mape = effectsize_cbe(p0_e1=p0_ce1, p0_e2=p0_e3,  
                         eff_e1=OR_ce1,
                         effm_e1 = "or",
                         eff_e2=OR3,
                         effm_e2 = "or",
                         effm_ce = "or",
                         rho=0
)$`Effect CE`

rho=c(0,0.1,0.2,0.3,0.4)

ss_cbe <- mapply(samplesize_cbe, p0_e1=p0_ce1,
                 p0_e2=p0_e3,
                 eff_e1=OR_ce1,
                 effm_e1="or",
                 eff_e2=OR3,
                 effm_e2="or",
                 effm_ce = "or",
                 rho=c(0,0.1,0.2,0.3,0.4),alpha=alpha, beta=beta) 

######################################
# Plot

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
ggsave("example_s.pdf",width = 100, height = 100, units = "mm")

######################################
# Power considerations

# Scenarios
p0_e1 = p0_ce1
p0_e2 = p0_e3
OR1 = OR_ce1
OR2 = OR3
corr_ce1 = "NO"
ss = mapply(samplesize_cbe, p0_e1=p0_ce1,
            p0_e2=p0_e3,
            eff_e1=OR_ce1,
            effm_e1="or",
            eff_e2=OR3,
            effm_e2="or",
            effm_ce = "or",
            rho=0,alpha=alpha, beta=beta)

scenarios1 = expand_grid(p0_e1,p0_e2,OR1,OR2,corr_ce1,ss)
scenarios=scenarios1

# Probabilities treat group
scenarios$p1_e1 = (scenarios$OR1*scenarios$p0_e1/(1-scenarios$p0_e1))/(1+(scenarios$OR1*scenarios$p0_e1/(1-scenarios$p0_e1)))
scenarios$p1_e2 = (scenarios$OR2*scenarios$p0_e2/(1-scenarios$p0_e2))/(1+(scenarios$OR2*scenarios$p0_e2/(1-scenarios$p0_e2)))

# 
corr = c(0, 0.1, 0.2, 0.3, 0.4)
dataset = expand_grid(scenarios, corr)

# head(dataset)
# 
# Calculate CE Probabilities
dataset$p0_ce = mapply(prob_cbe, p_e1=dataset$p0_e1, p_e2=dataset$p0_e2, rho=dataset$corr)
dataset$p1_ce = mapply(prob_cbe, p_e1=dataset$p1_e1, p_e2=dataset$p1_e2, rho=dataset$corr) 
dataset$OR_ce = (dataset$p1_ce/(1-dataset$p1_ce))/(dataset$p0_ce/(1-dataset$p0_ce))

head(dataset)

dataset$Power<-c()
dataset$Power_CE<-c()
dataset$Power_RE<-c()
dataset$decision<-c()

nsim=100000
# type i and ii errors
z.alpha <- qnorm(1-alpha,0,1)  
z.beta <-  qnorm(1-beta,0,1) 

for(i in 1:dim(dataset)[1]){
  aux <- rowSums(replicate(nsim,eselectsim(ss_arm=round(dataset$ss[i]/2),
                                           p0_e1=dataset$p0_e1[i],OR1=dataset$OR1[i],
                                           p0_e2=dataset$p0_e2[i],OR2=dataset$OR2[i],
                                           p0_ce=dataset$p0_ce[i],p_init=0.5,
                                           criteria="SS",H0_e1=FALSE,H0_e2=FALSE,SS_r=TRUE))< c(-z.alpha,1))/nsim 
  
  
  dataset$Power[i]<- aux[1]
  dataset$decision[i]<- 1-aux[2]
  print(i)
}

for(i in 1:dim(dataset)[1]){
  dataset$Power_CE[i] <- sum(replicate(nsim,f_OR(dataset$ss[i]/2,dataset$p0_ce[i],dataset$p1_ce[i]))< - z.alpha)/nsim
  dataset$Power_RE[i] <- sum(replicate(nsim,f_OR(dataset$ss[i]/2,dataset$p0_e1[i],dataset$p1_e1[i]))< - z.alpha)/nsim
  print(i)
}

summary(dataset)


dataset$pe <- ifelse(dataset$decision>1,"RE","CE")
# 
power_data <- data.frame(Power=c(dataset$Power,dataset$Power_CE,dataset$Power_RE),
                         Design=c(rep("AD",length(dataset$Power)), 
                                  rep("CD",length(dataset$Power_CE)),
                                  rep("RD",length(dataset$Power_RE))),
                         Correlation=c(dataset$corr_ce1,dataset$corr_ce1,dataset$corr_ce1),
                         corr=c(dataset$corr,dataset$corr,dataset$corr) 
)

power_data$Design <- factor(power_data$Design, levels=c("CD","RD","AD"))
power_data$Correlation <- as.factor(power_data$Correlation)
# 
# # 
ggplot(power_data,aes(x = corr, y = Power, color = Design#, linetype = Correlation
)) + geom_point() +geom_line(size=0.8) + xlab("Correlation") + ylim(c(0.7,0.85)) + scale_color_manual(breaks = c("CD","RD","AD"), values=c("#D9717D", "#4DB6D0", "#BECA55")) 
ggsave("example_power_s.pdf",width = 100, height = 100, units = "mm") 





