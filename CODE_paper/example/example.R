
##################################################################################
# Research project Vienna - Bcn
# Endpoint selection and sample size reassessment for composite binary endpoints
# Simulation study 
##################################################################################

rm(list = ls())

library(CompAREdesign)
library(tidyverse)
library(gridExtra)

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

OR1=0.66
OR2=0.54
OR3=0.64

alpha=0.05
beta=0.2

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




######################################
# Relevant endpoint
p0_ce1 = prob_cbe(p_e1=p0_e2, p_e2=p0_e3, rho=0)
p0_ce1_corr = prob_cbe(p_e1=p0_e2, p_e2=p0_e3, rho=0.5)

OR_ce1 = effectsize_cbe(p0_e1=p0_e2, p0_e2=p0_e3,  
  eff_e1=OR2,
  effm_e1 = "or",
  eff_e2=OR3,
  effm_e2 = "or",
  effm_ce = "or",
  rho=0
)$`Effect CE`

OR_ce1_corr = effectsize_cbe(p0_e1=p0_e2, p0_e2=p0_e3,  
                        eff_e1=OR2,
                        effm_e1 = "or",
                        eff_e2=OR3,
                        effm_e2 = "or",
                        effm_ce = "or",
                        rho=0.5
)$`Effect CE`

n1 = samplesize_OR(p0=p0_ce1, OR=OR_ce1, alpha=alpha, beta=beta)
n1_corr = samplesize_OR(p0=p0_ce1_corr, OR=OR_ce1_corr, alpha=alpha, beta=beta)
  
# Composite endpoint
# corr ounds
lower_corr(p_e1=p0_ce1, p_e2=p0_e3)
upper_corr(p_e1=p0_ce1, p_e2=p0_e3)

lower_corr(p_e1=p0_ce1_corr, p_e2=p0_e3)
upper_corr(p_e1=p0_ce1_corr, p_e2=p0_e3)


prob_cbe(p_e1=p0_ce1, p_e2=p0_e3, rho=0)
prob_cbe(p_e1=p0_ce1_corr, p_e2=p0_e3, rho=0)

rho=c(0,0.1,0.2,0.3,0.4)

ss_cbe <- mapply(samplesize_cbe, p0_e1=p0_ce1,
       p0_e2=p0_e3,
       eff_e1=OR_ce1,
       effm_e1="or",
       eff_e2=OR3,
       effm_e2="or",
       effm_ce = "or",
       rho=c(0,0.1,0.2,0.3,0.4),alpha=alpha, beta=beta) 

ss_cbe_corr <- mapply(samplesize_cbe, p0_e1=p0_ce1_corr,
                 p0_e2=p0_e3,
                 eff_e1=OR_ce1_corr,
                 effm_e1="or",
                 eff_e2=OR3,
                 effm_e2="or",
                 effm_ce = "or",
                 rho=c(0,0.1,0.2,0.3,0.4),alpha=alpha, beta=beta) 
# 
data_plot = data.frame(ss=c(ss_cbe,
                            ss_cbe_corr,
                            rep(n1,length(ss_cbe)),
                            rep(n1_corr,length(ss_cbe))),
                       rho=rep(rho,4),
                       Endpoint=c(rep("CE",2*length(ss_cbe)),rep("RE",2*length(ss_cbe))
                             ),
                       Correlation=c(rep("NO",length(ss_cbe)),rep("YES",length(ss_cbe)),
                              rep("NO",length(ss_cbe)),rep("YES",length(ss_cbe))
                              )
                       )
data_plot$Endpoint = as.factor(data_plot$Endpoint)
data_plot$Correlation = as.factor(data_plot$Correlation)

ggplot(data_plot,aes(x = rho, y = ss, color=Endpoint, linetype = Correlation)) + geom_point() +geom_line(size=0.8) + xlab("Correlation") + ylab("Sample size") 
ggsave("example_ss.pdf",width = 100, height = 100, units = "mm")

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
# 
p0_e1 = p0_ce1_corr
p0_e2 = p0_e3
OR1 = OR_ce1_corr
OR2 = OR3
corr_ce1 = "YES"
ss = mapply(samplesize_cbe, p0_e1=p0_ce1_corr,
            p0_e2=p0_e3,
            eff_e1=OR_ce1_corr,
            effm_e1="or",
            eff_e2=OR3,
            effm_e2="or",
            effm_ce = "or",
            rho=0,alpha=alpha, beta=beta)

scenarios2 = expand_grid(p0_e1,p0_e2,OR1,OR2,corr_ce1,ss)
scenarios=rbind(scenarios1,scenarios2) 

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
# dataset$Power_RE<-c()
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
  # dataset$Test_Power_RE[i] <- sum(replicate(nsim,f_OR(dataset$ss[i]/2,dataset$p0_e1[i],dataset$p0_e1[i]))< - z.alpha)/nsim
  print(i)
}

summary(dataset)


dataset$pe <- ifelse(dataset$decision>1,"RE","CE")
# 
power_data <- data.frame(Power=c(dataset$Power,dataset$Power_CE),
                         Design=c(rep("AD",length(dataset$Power)), 
                                  rep("CE",length(dataset$Power_CE))),
                         Correlation=c(dataset$corr_ce1,dataset$corr_ce1),
                         corr=c(dataset$corr,dataset$corr) 
)

power_data$Design <- as.factor(power_data$Design)
power_data$Correlation <- as.factor(power_data$Correlation)
# 
# # 
ggplot(power_data,aes(x = corr, y = Power, color = Design, linetype = Correlation)) + geom_point() +geom_line(size=0.8) + xlab("Correlation") + ylim(c(0.7,0.85))
# + labs(color = "Decision\n", linetype = "Correlation\n")
# # + ylab("Sample size")
ggsave("example_power.pdf",width = 100, height = 100, units = "mm")
