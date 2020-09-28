
##################################################################################
# Research project Vienna - Bcn
# Endpoint selection and sample size reassessment for composite binary endpoints
# Simulation study 
##################################################################################

rm(list = ls())

# setwd("C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/CBE_selection")
# source('C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/CBE_selection/Functions.R') 
# source('C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/CBE_selection/sim_functions_ce.R')

setwd("C:/Users/Marta/Nextcloud/Gitkraken/CBE_selection")
source('C:/Users/Marta/Nextcloud/GitKraken/CBE_selection/Functions.R') 
source('C:/Users/Marta/Nextcloud/GitKraken/CBE_selection/sim_functions_ce.R') 

##################################################################################
# Differences with respect to previous versions
##################################################################################

# - consider less scenarios for OR2, p1, p2
# - add scenarios correlation between components
# - computation sample size according to E1

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
p0_e2 = c(0.1, 0.25)
OR1 = c(0.6, 0.7, 0.8)
OR2 = c(0.75, 0.8)

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

corr = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)

dataset = expand_grid(scenarios, corr) 
dataset = subset(dataset, dataset$corr < dataset$max_corr & dataset$corr > dataset$min_corr)

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

# clean
rm(OR1,OR2,p0_e1,p0_e2,corr,scenarios)

# Vector empirical powers and significance level
dataset$Test_Power_CE <- NA
dataset$Test_Power_RE <- NA
dataset$Test_Power_ES <- NA
dataset$Test_Power_ES_ub <- NA
dataset$decision_ES <- NA 
dataset$decision_ES_ub <- NA 

dataset$Test_Reject_CE <- NA
dataset$Test_Reject_RE <- NA
dataset$Test_Reject_ES <- NA
dataset$Test_Reject_ES_ub <- NA

#########################################
# Simulations
#########################################

set.seed(4123)

# nsim: number of simulations
nsim = 100000 

# type i and ii errors
z.alpha <- qnorm(1-alpha,0,1)  
z.beta <-  qnorm(1-beta,0,1) 


t0=Sys.time() 

#########################################

for(i in 1:dim(dataset)[1]){
  dataset$Test_Power_CE[i] <- sum(replicate(nsim,f_OR(dataset$samplesize_e1[i]/2,dataset$p0_ce[i],dataset$p1_ce[i]))< - z.alpha)/nsim
  dataset$Test_Power_RE[i] <- sum(replicate(nsim,f_OR(dataset$samplesize_e1[i]/2,dataset$p0_e1[i],dataset$p1_e1[i]))< - z.alpha)/nsim
  print(i)
} 

#  

for(i in 1:dim(dataset)[1]){
  aux <- rowSums(replicate(nsim,eselection_b(samplesize=dataset$samplesize_e1[i]/2,
                                   p0_e1=dataset$p0_e1[i],p1_e1=dataset$p1_e1[i],
                                   OR1=dataset$OR1[i],
                                   p0_e2=dataset$p0_e2[i],p1_e2=dataset$p1_e2[i],
                                   OR2=dataset$OR2[i],
                                   p0_ce=dataset$p0_ce[i],p1_ce=dataset$p1_ce[i]))< c(-z.alpha,1))/nsim
  
                                                       
                                         
  dataset$Test_Power_ES[i]<- aux[1]
  dataset$decision_ES[i]<- 1-aux[2]
  print(i)
}

#  

for(i in 1:dim(dataset)[1]){ 
  aux <- rowSums(replicate(nsim,eselection_ub(samplesize=dataset$samplesize_e1[i]/2,
                                         p0_e1=dataset$p0_e1[i],p1_e1=dataset$p1_e1[i],
                                         OR1=dataset$OR1[i],
                                         p0_e2=dataset$p0_e2[i],p1_e2=dataset$p1_e2[i],
                                         OR2=dataset$OR2[i],
                                         p0_ce=dataset$p0_ce[i],p1_ce=dataset$p1_ce[i]))< c(-z.alpha,1))/nsim
                  
  dataset$Test_Power_ES_ub[i]<- aux[1]
  dataset$decision_ES_ub[i]<- 1-aux[2]
  print(i)
} 
save.image("C:/Users/Marta/Dropbox/C5/Scripts/GitKraken/CBE_selection/results/results_H0_False.RData") 
#########################################


set.seed(2314)

for(i in 1:dim(dataset)[1]){
  dataset$Test_Reject_CE[i] <- sum(replicate(nsim,f_OR(dataset$samplesize_e1[i]/2,dataset$p0_ce[i],dataset$p0_ce[i]))< - z.alpha)/nsim
  dataset$Test_Reject_RE[i] <- sum(replicate(nsim,f_OR(dataset$samplesize_e1[i]/2,dataset$p0_e1[i],dataset$p0_e1[i]))< - z.alpha)/nsim
  print(i)
}


#########################################  

for(i in 1:dim(dataset)[1]){
  dataset$Test_Reject_ES[i] <- sum(replicate(nsim,
                                             eselection_b(samplesize=dataset$samplesize_e1[i]/2,
                                                   p0_e1=dataset$p0_e1[i],p1_e1=dataset$p0_e1[i],
                                                   OR1=dataset$OR1[i],
                                                   p0_e2=dataset$p0_e2[i],p1_e2=dataset$p0_e2[i],
                                                   OR2=dataset$OR2[i],
                                                   p0_ce=dataset$p0_ce[i],p1_ce=dataset$p0_ce[i])[1])< - z.alpha)/nsim
  # tryCatch(, error=function(e){NA})
  print(i)
}

# 


for(i in 1:dim(dataset)[1]){
  dataset$Test_Reject_ES_ub[i] <- sum(replicate(nsim,
                                                eselection_ub(samplesize=dataset$samplesize_e1[i]/2,
                                                         p0_e1=dataset$p0_e1[i],p1_e1=dataset$p0_e1[i],
                                                         OR1=dataset$OR1[i],
                                                         p0_e2=dataset$p0_e2[i],p1_e2=dataset$p0_e2[i],
                                                         OR2=dataset$OR2[i],
                                                         p0_ce=dataset$p0_ce[i],p1_ce=dataset$p0_ce[i])[1])< - z.alpha)/nsim
  print(i)
}

t1=Sys.time()-t0  

save.image("C:/Users/Marta/Dropbox/C5/Scripts/GitKraken/CBE_selection/results/results_H0_True.RData")
# save.image("C:/Users/Marta/Dropbox/C5/Scripts/GitKraken/CBE_selection/results/results_unblinded.RData")



#########################################  

######################################### 
#########################################  
# results
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(tidyverse) 
load("C:/Users/Marta/Nextcloud/Gitkraken/CBE_selection/scenarios.RData")

set.seed(42)

# nsim: number of simulations
nsim = 100

t0=Sys.time() 

corr_est_ub <- data.frame(numeric(nsim+1))
corr_est_b <- data.frame(numeric(nsim+1))

v_ub <- c()
v_b <- c()

for(i in 1:dim(dataset)[1]){
  v_ub <- replicate(nsim, estimation_ub(samplesize=dataset$samplesize_e1[i]/2,
                                                  p0_e1=dataset$p0_e1[i],p1_e1=dataset$p1_e1[i],
                                                  OR1=dataset$OR1[i],
                                                  p0_e2=dataset$p0_e2[i],p1_e2=dataset$p1_e2[i],
                                                  OR2=dataset$OR2[i],
                                                  p0_ce=dataset$p0_ce[i],p1_ce=dataset$p1_ce[i])[1])
  
  corr_est_ub[,i] <- c(dataset$corr[i],v_ub)
  
  v_b <- replicate(nsim, estimation_b(samplesize=dataset$samplesize_e1[i]/2,
                                      p0_e1=dataset$p0_e1[i],p1_e1=dataset$p1_e1[i],
                                      OR1=dataset$OR1[i],
                                      p0_e2=dataset$p0_e2[i],p1_e2=dataset$p1_e2[i],
                                      OR2=dataset$OR2[i],
                                      p0_ce=dataset$p0_ce[i],p1_ce=dataset$p1_ce[i])[1])
  
  corr_est_b[,i] <- c(dataset$corr[i],v_b)
  
  print(i)
}

t1=Sys.time()-t0  
#  
save.image("C:/Users/Marta/Nextcloud/Gitkraken/CBE_selection/results/corr_estimation.RData")

corr_est_b <- corr_est_b %>%
              tibble::rownames_to_column() %>%  
              pivot_longer(-rowname) %>% 
              pivot_wider(names_from=rowname, values_from=value) 
corr_est_b <- corr_est_b[-1]

corr_est_ub <- corr_est_ub %>%
                tibble::rownames_to_column() %>%  
                pivot_longer(-rowname) %>% 
                pivot_wider(names_from=rowname, values_from=value) 
              
corr_est_ub <- corr_est_ub[-1]

# 

windows()
p <- list() 
enum = 1 
it = 1
i =1
for(i in 1:max(dataset$scenario)){
  
  sub=subset(dataset,dataset$scenario==i)
  sub$corr=as.factor(sub$corr)
  
  dataggplot <- data.frame(corr=as.factor(corr_est_b[1:9,1]), corr_est=corr_est_b[2:9,])
  
  p[[enum]] <-ggplot(dataggplot,aes(x=corr, y=corr_est))  + geom_boxplot() + ggtitle(paste("Scenario", dataset$scenario[it], "\n (p1,p2,OR1,OR2) \n=(", dataset$p0_e1[it],",",dataset$p0_e2[it],",",dataset$OR1[it],",",dataset$OR2[it],")"))
  
  p[[enum+1]] <-ggplot(sub, aes(x=corr, y=diff_corr_b))  + geom_boxplot() + ggtitle(paste("Scenario", dataset$scenario[it], "\n (p1,p2,OR1,OR2) \n=(", dataset$p0_e1[it],",",dataset$p0_e2[it],",",dataset$OR1[it],",",dataset$OR2[it],")")) 
  
  # stable <- data.frame(Corr=sub$corr,ARE=round(sub$ARE,2),"%CE"=round(100*sub$decision_ES,2),"%CE_u"=round(100*sub$decision_ES_ub,2),check.names=FALSE)  
  # 
  # p[[enum+4]] <- ggtexttable(stable, rows = NULL, theme = ttheme(base_style = "default", base_size = 9))
  
  enum=enum+3 
  it <- it + dim(subset(dataset,dataset$scenario==i))
}
marrangeGrob(p,ncol=2,nrow=1,top=NULL) 








