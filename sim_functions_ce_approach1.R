

# 2x2 Table of probabilities
s1_group0 = phi_group0
s2_group0 = p1_group0 - phi_group0
s3_group0 = p2_group0 - phi_group0
s4_group0 = 1 - s1_group0 - s2_group0 - s3_group0

  

##############################################################  
i=1
p0_e1 = dataset$p0_e1[i]
p1_e1 = dataset$p1_e1[i]
p0_e2 = dataset$p0_e2[i]
p1_e2 = dataset$p1_e2[i]
p0_ce = dataset$p0_ce[i]
p1_ce = dataset$p1_ce[i]
samplesize=500
# dataset$corr[i]



s1_group0 = p0_e1+p0_e2-p0_ce
s2_group0 = p0_e1 - (p0_e1+p0_e2-p0_ce)
s3_group0 = p0_e2 - (p0_e1+p0_e2-p0_ce)
s4_group0 = 1 - s1_group0 - s2_group0 - s3_group0

sm = rmultinom(1,500,c(s1_group0,s2_group0,s3_group0,s4_group0))
corr_est =(sm[1]*sm[4]-(sm[2]*sm[3]))/sqrt((sm[1]+sm[3])*(sm[2]+sm[4])*(sm[1]+sm[2])*(sm[3]+sm[4]))
(corr_est)

##############################################################  

f_ES <- function(samplesize,p0_e1,p1_e1,p0_e2,p1_e2,p0_ce,p1_ce){
  
  rmultinom(2,500,c(0.75,0.25))
  
  # CONTROL GROUP
  # Endpoint 1
  R = runif(samplesize)
  X1_0 = sum(R<p0_e1 & R>0)
  # phat_group0 = X1_0/samplesize
  
  # Endpoint 2
  R = runif(X1_0)
  p_aux=(p0_e1+p0_e2-p0_ce)/p0_e1
  X21_0 = sum(R<p_aux & R>0)
  # phat_group0 = X21_0/samplesize
  
  R = runif(samplesize-X1_0)
  p_aux=(p0_ce-p0_e1)/(1-p0_e1)
  X20_0 = sum(R<p_aux & R>0)
  # phat_group0 = X20_0/samplesize
  
  # INTERVENTION GROUP
  # Endpoint 1
  R = runif(samplesize)
  X1_1 = sum(R<p1_e1 & R>0)
  # phat_group0 = X1_1/samplesize
  
  # Endpoint 2
  R = runif(X1_1)
  p_aux=(p1_e1+p1_e2-p1_ce)/p1_e1
  X21_1 = sum(R<p_aux & R>0)
  # phat_group0 = X21_1/samplesize
  
  R = runif(samplesize-X1_1)
  p_aux=(p1_ce-p1_e1)/(1-p1_e1)
  X20_1 = sum(R<p_aux & R>0)
  # phat_group0 = X20_1/samplesize
  
  X1 = c(rep(1,X1_0),rep(0,samplesize-X1_0),rep(1,X1_1),rep(0,samplesize-X1_1))
  X2 = c(rep(1,X21_0),rep(1,X20_0),rep(0,samplesize-X21_0-X20_0),rep(1,X21_1),rep(1,X20_1),rep(0,samplesize-X21_1-X20_1))
  
  # cor(X1,X2)
  
  sum(X1)/samplesize
  sum(X2)/samplesize
  
  
  
  # X1 = c(rep(1,X1_0),rep(0,samplesize-X1_0))
  # X2 = c(rep(1,X21_0),rep(1,X20_0),rep(0,samplesize-X21_0-X20_0))
  
  
  
  # test odds ratio
  TestOR_unpooled = log((phat_group1/(1-phat_group1))/(phat_group0/(1-phat_group0)))*((1/(phat_group0*(1-phat_group0))+ 1/(phat_group1*(1-phat_group1)))/samplesize)^(-1/2)
  
  return(TestOR_unpooled)
}