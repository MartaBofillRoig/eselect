#####################################################################################
# Sample size calculations in terms of odds ratio  

samplesize_OR <- function(p0, OR, alpha=0.05, beta=0.2, Unpooled="Unpooled Variance"){
  z.alpha <- qnorm(1-alpha,0,1)  
  z.beta <-  qnorm(1-beta,0,1)
  
  p1 = (OR*p0/(1-p0))/(1+(OR*p0/(1-p0)))
  
  if(Unpooled=="Unpooled Variance"){
    # sample size per group
    n1 = ((z.alpha+z.beta)/(log(OR)))^2*( 1/(p0*(1-p0)) + 1/(p1*(1-p1)) )  
  }else{
    p = (p1 + p0)/2
    # sample size per group
    n1 = ((z.alpha* sqrt(2/(p*(1-p))) + z.beta* sqrt(1/(p0*(1-p0)) + 1/(p1*(1-p1))))/(log(OR)))^2 
  }
  
  n = 2*n1  
  
  return(n)
}


############################################################## 
# f_OR; simulates binary outcomes and computes the statistic

f_OR <- function(samplesize,p0,p1){
  
  # estimation of the probability of observing the composite event in control group
  R = runif(samplesize)
  phat_group0 = 1 - sum(R>=p0 & R<1)/samplesize
  
  # estimation of the probability of observing the composite event in test group
  R = runif(samplesize)
  phat_group1 = 1 - sum(R>=p1 & R<1)/samplesize
  
  # test odds ratio 
  TestOR_unpooled = test_f(OR=OR_function(p0=phat_group0,p1=phat_group1),p0=phat_group0,n=samplesize)
  
  return(TestOR_unpooled)
} 


############################################################## 
# test_f: computes the statistical test (OR)

test_f <- function(OR,p0,n){
  
  # OR: odds ratio
  # p0: probability under control group
  # n: sample size per group
  
  p1 <- (OR*p0/(1-p0))/(1+(OR*p0/(1-p0)))
  den <- 1/(n*(p0*(1-p0))) + 1/(n*(p1*(1-p1)))
  
  test <- (log(OR))/sqrt(den)
  return(test)
  
} 

############################################################## 
# test_me: computes the statistical tests (OR) for two endpoints

test_me <- function(OR1,p0_e1,OR2,p0_e2,n){
  
  # OR: (estimated) odds ratio 
  # p0: (estimated) probability under control group
  # n: sample size per group
  
  p1 <- (OR1*p0_e1/(1-p0_e1))/(1+(OR1*p0_e1/(1-p0_e1)))
  den <- 1/(n*(p0_e1*(1-p0_e1))) + 1/(n*(p1*(1-p1)))
  test1 <- (log(OR1))/sqrt(den)
  
  p1 <- (OR2*p0_e2/(1-p0_e2))/(1+(OR2*p0_e2/(1-p0_e2)))
  den <- 1/(n*(p0_e2*(1-p0_e2))) + 1/(n*(p1*(1-p1)))
  test2 <- (log(OR2))/sqrt(den)
  
  return(list(Test1=test1,Test2=test2))
  
} 

############################################################## 
# simtest_me: 
# simulates two correlated binary endpoints and
# computes the statistical test (OR) for two endpoints

simtest_me <- function(OR1,p0_e1,OR2,p0_e2,n){ 
  
  # control group
  sm0 = f_sim(samplesize=ss_arm,p_e1=p0_e1,p_e2=p0_e2,p_ce=p0_ce)
  # intervention group
  sm1 = f_sim(samplesize=ss_arm,p_e1=p1_e1,p_e2=p1_e2,p_ce=p1_ce)
  
  # estimated probabilities
  phat0_e1 = (sm0[1]+sm0[2])/samplesize
  phat0_e2 = (sm0[1]+sm0[3])/samplesize
  
  phat1_e1 = (sm1[1]+sm1[2])/samplesize
  phat1_e2 = (sm1[1]+sm1[3])/samplesize
  
  test_univ = test_fme(OR1=OR_function(p0=phat0_e1, p1=phat1_e1),p0_e1=phat0_e1,
                       OR2=OR_function(p0=phat0_e2, p1=phat1_e2),p0_e2=phat0_e2,
                       n=ss_arm)
  
  return(list(Test1=test1,Test2=test2))
  
} 


##############################################################  
# OR_function: computes the odds ratio

OR_function<- function(p0, p1){
  
  OR<- (p1/(1-p1))/(p0/(1-p0))
  return(OR)
}

##############################################################  
# fun_p0: computes the probability under the control group
# based on the pooled probability and the odds ratio

# p: pooled probability
# l: odds ratio

fun_p0 <- function(p,l){
  
  p = p*2
  
  sol1= -(1 + l + p - l*p + sqrt(4*(-1 + l)*p + (1 + l + p - l*p)^2))/(2*(-1 + l))
  sol2= -(1 + l + p - l*p - sqrt(4*(-1 + l)*p + (1 + l + p - l*p)^2))/(2*(-1 + l))
  
  if(sol2>0 & sol2<1){
    p0=sol2
  }else{
    p0=sol1
  }
  
  return(p0)  
}
# EXAMPLE
# p=p1_e1+p0_e1
# l=OR1
# fun_p0(p=p,l=l)

##############################################################  
# corr_rest: computes the correlation restrictions

# phat0_e1,phat0_e2,phat1_e1,phat1_e2: estimated probabilities
# p0_e1,p0_e2,p1_e1,p1_e2: assumed probabilities (design)
# OR1,OR2: expected effect sizes

corr_rest_ub <- function(phat0_e1,phat0_e2,phat1_e1,phat1_e2,p0_e1,p0_e2,p1_e1,p1_e2,OR1,OR2){ 
  # correlation restrictions 
  update_uppcorr0=upper_corr(phat0_e1,phat0_e2)
  update_uppcorr1=upper_corr(phat1_e1,phat1_e2) 
  uppcorr0=upper_corr(p0_e1,p0_e2)
  uppcorr1=upper_corr(p1_e1,p1_e2)
  uppcorr12=upper_corr((OR1*phat0_e1/(1-phat0_e1))/(1+(OR1*phat0_e1/(1-phat0_e1))),
                       (OR2*phat0_e2/(1-phat0_e2))/(1+(OR2*phat0_e2/(1-phat0_e2))))
  
  upp = min(update_uppcorr0,update_uppcorr1,uppcorr0,uppcorr1,uppcorr12)
  
  update_lowcorr0= lower_corr(phat0_e1,phat0_e2)  
  update_lowcorr1=lower_corr(phat1_e1,phat1_e2)
  lowcorr0= lower_corr(p0_e1,p0_e2)  
  lowcorr1=lower_corr(p1_e1,p1_e2)
  lowcorr12=lower_corr((OR1*phat0_e1/(1-phat0_e1))/(1+(OR1*phat0_e1/(1-phat0_e1))),
                       (OR2*phat0_e2/(1-phat0_e2))/(1+(OR2*phat0_e2/(1-phat0_e2))))
  
  low = max(update_lowcorr0,update_lowcorr1,lowcorr0,lowcorr1,lowcorr12) 
  
  return(c(low,upp))  
}

##############################################################  
# corr_rest: computes the correlation restrictions

# phat0_e1,phat0_e2,phat1_e1,phat1_e2: estimated probabilities
# p0_e1,p0_e2,p1_e1,p1_e2: assumed probabilities (design)
# OR1,OR2: expected effect sizes

corr_rest_b <- function(phat0_e1,phat0_e2,p0_e1,p0_e2,p1_e1,p1_e2,OR1,OR2){ 
  # correlation restrictions  
  update_uppcorr0=upper_corr(phat0_e1,phat0_e2) 
  uppcorr0=upper_corr(p0_e1,p0_e2)
  uppcorr1=upper_corr(p1_e1,p1_e2)
  uppcorr12=upper_corr((OR1*phat0_e1/(1-phat0_e1))/(1+(OR1*phat0_e1/(1-phat0_e1))),
                       (OR2*phat0_e2/(1-phat0_e2))/(1+(OR2*phat0_e2/(1-phat0_e2))))
  
  upp = min(update_uppcorr0,uppcorr0,uppcorr1,uppcorr12)
  
  update_lowcorr0= lower_corr(phat0_e1,phat0_e2)  
  lowcorr0= lower_corr(p0_e1,p0_e2)  
  lowcorr1=lower_corr(p1_e1,p1_e2)
  lowcorr12=lower_corr((OR1*phat0_e1/(1-phat0_e1))/(1+(OR1*phat0_e1/(1-phat0_e1))),
                       (OR2*phat0_e2/(1-phat0_e2))/(1+(OR2*phat0_e2/(1-phat0_e2))))
  
  low = max(update_lowcorr0,lowcorr0,lowcorr1,lowcorr12)  
  
  return(c(low,upp))  
}




##############################################################  
# f_sim:
# simulation binary data two outcomes 

# returns 2x2 table
# s1+s2: num X1
# s1+s3: num X2

f_sim <- function(samplesize,p_e1,p_e2,p_ce){ 
  
  # 2x2 table
  s1_group = p_e1 + p_e2 - p_ce
  s2_group = ifelse(p_ce-p_e2>0, p_ce-p_e2, 0)#p_ce-p_e2  
  s3_group = ifelse(p_ce-p_e1>0, p_ce-p_e1, 0)#p_ce-p_e1  
  s4_group = 1- p_ce  
  
  data = rmultinom(n=1,size=round(samplesize),prob=c(s1_group,s2_group,s3_group,s4_group)) 
  
  return(data)
}


##############################################################  
# estimation_ub:
# simulation
# estimation correlation unblinded without previous info corr/CE 

estimation_ub <- function(samplesize,p0_e1,p1_e1,OR1,p0_e2,p1_e2,OR2,p0_ce,p1_ce){  
  
  total_ss = samplesize*2
  
  # control group
  sm0 = f_sim(samplesize=samplesize,p_e1=p0_e1,p_e2=p0_e2,p_ce=p0_ce)
  
  # intervention group
  sm1 = f_sim(samplesize=samplesize,p_e1=p1_e1,p_e2=p1_e2,p_ce=p1_ce)
  
  # pooled sample
  sm = sm0 + sm1
  
  # estimated probabilities 
  phat0_e1 = (sm0[1]+sm0[2])/samplesize
  phat0_e2 = (sm0[1]+sm0[3])/samplesize
  phat0_ce = 1-(sm0[4])/samplesize
  
  phat1_e1 = (sm1[1]+sm1[2])/samplesize
  phat1_e2 = (sm1[1]+sm1[3])/samplesize
  phat1_ce = 1-(sm1[4])/samplesize
  
  # estimated correlation
  corrhat0 = ((phat0_e1+phat0_e2-phat0_ce)-phat0_e1*phat0_e2)/sqrt(phat0_e1*(1-phat0_e1)*phat0_e2*(1-phat0_e2))
  corrhat1 = ((phat1_e1+phat1_e2-phat1_ce)-phat1_e1*phat1_e2)/sqrt(phat1_e1*(1-phat1_e1)*phat1_e2*(1-phat1_e2))
  
  corrhat = (corrhat0+corrhat1)/2
  
  # correlation restrictions 
  update_uppcorr0=upper_corr(phat0_e1,phat0_e2)
  update_uppcorr1=upper_corr(phat1_e1,phat1_e2) 
  uppcorr0=upper_corr(p0_e1,p0_e2)
  uppcorr1=upper_corr(p1_e1,p1_e2)
  uppcorr12=upper_corr((OR1*phat0_e1/(1-phat0_e1))/(1+(OR1*phat0_e1/(1-phat0_e1))),
                       (OR2*phat0_e2/(1-phat0_e2))/(1+(OR2*phat0_e2/(1-phat0_e2))))
  
  upp = min(update_uppcorr0,update_uppcorr1,uppcorr0,uppcorr1,uppcorr12)
  
  update_lowcorr0= lower_corr(phat0_e1,phat0_e2)  
  update_lowcorr1=lower_corr(phat1_e1,phat1_e2)
  lowcorr0= lower_corr(p0_e1,p0_e2)  
  lowcorr1=lower_corr(p1_e1,p1_e2)
  lowcorr12=lower_corr((OR1*phat0_e1/(1-phat0_e1))/(1+(OR1*phat0_e1/(1-phat0_e1))),
                       (OR2*phat0_e2/(1-phat0_e2))/(1+(OR2*phat0_e2/(1-phat0_e2))))
  
  low = max(update_lowcorr0,update_lowcorr1,lowcorr0,lowcorr1,lowcorr12)
  
  corrhat_c = corrhat
  if(corrhat > upp){
    corrhat_c = upp
  }
  if(corrhat < low){
    corrhat_c = low
  } 
  
  return(list=c(corrhat=corrhat,corrhat_c=corrhat_c,phat0_e1=phat0_e1,phat0_e2=phat0_e2))
}  


##############################################################  
# estimation_b:
# simulation
# estimation correlation blinded without previous info corr/CE 

estimation_b <- function(samplesize,p0_e1,p1_e1,OR1,p0_e2,p1_e2,OR2,p0_ce,p1_ce){  
  total_ss = samplesize*2
  # control group
  s1_group0 = p0_e1 + p0_e2 - p0_ce
  s2_group0 = p0_ce-p0_e2  
  s3_group0 = p0_ce-p0_e1  
  s4_group0 = 1- p0_ce  
  
  sm0 = rmultinom(1,samplesize,c(s1_group0,s2_group0,s3_group0,s4_group0)) 
  
  # intervention group
  s1_group1 = p1_e1+p1_e2-p1_ce
  s2_group1 = ifelse(p1_ce-p1_e2>0, p1_ce-p1_e2, 0)  
  s3_group1 = ifelse(p1_ce-p1_e1>0, p1_ce-p1_e1, 0)   
  s4_group1 = 1- p1_ce  
  
  sm1 = rmultinom(1,samplesize,c(s1_group1,s2_group1,s3_group1,s4_group1))
  
  # pooled sample
  sm = sm0 + sm1
  
  # estimated probabilities
  phat_e1 = (sm[1]+sm[2])/total_ss
  phat_e2 = (sm[1]+sm[3])/total_ss
  phat_ce = 1-(sm[4])/total_ss
  
  # 
  phat0_e1 = fun_p0(p=phat_e1,l=OR1)
  phat0_e2 = fun_p0(p=phat_e2,l=OR2) 
  
  # 
  phat1_e1 = (OR1*phat0_e1/(1-phat0_e1))/(1+(OR1*phat0_e1/(1-phat0_e1)))
  phat1_e2 = (OR2*phat0_e2/(1-phat0_e2))/(1+(OR2*phat0_e2/(1-phat0_e2)))
  
  # estimated correlation
  corrhat =  (phat_ce - (samplesize/total_ss)*(1 - (1-phat0_e1)*(1-phat0_e2)) - (samplesize/total_ss)*(1-(1-phat1_e1)*(1-phat1_e2)))/(-(samplesize/total_ss)*sqrt(phat0_e1*phat0_e2*(1-phat0_e1)*(1-phat0_e2))-(samplesize/total_ss)*sqrt(phat1_e1*phat1_e2*(1-phat1_e1)*(1-phat1_e2)))
  
  # correlation restrictions  
  update_uppcorr0=upper_corr(phat0_e1,phat0_e2)
  # update_uppcorr1=upper_corr(phat1_e1,phat1_e2) 
  uppcorr0=upper_corr(p0_e1,p0_e2)
  uppcorr1=upper_corr(p1_e1,p1_e2)
  uppcorr12=upper_corr((OR1*phat0_e1/(1-phat0_e1))/(1+(OR1*phat0_e1/(1-phat0_e1))),
                       (OR2*phat0_e2/(1-phat0_e2))/(1+(OR2*phat0_e2/(1-phat0_e2))))
  
  upp = min(update_uppcorr0,uppcorr0,uppcorr1,uppcorr12)
  
  update_lowcorr0= lower_corr(phat0_e1,phat0_e2)  
  # update_lowcorr1=lower_corr(phat1_e1,phat1_e2)
  lowcorr0= lower_corr(p0_e1,p0_e2)  
  lowcorr1=lower_corr(p1_e1,p1_e2)
  lowcorr12=lower_corr((OR1*phat0_e1/(1-phat0_e1))/(1+(OR1*phat0_e1/(1-phat0_e1))),
                       (OR2*phat0_e2/(1-phat0_e2))/(1+(OR2*phat0_e2/(1-phat0_e2))))
  
  low = max(update_lowcorr0,lowcorr0,lowcorr1,lowcorr12) 
  
  corrhat_c = corrhat
  if(corrhat > upp){
    corrhat_c = upp
  }
  if(corrhat < low){
    corrhat_c = low
  } 
  
  return(list=c(corrhat=corrhat,corrhat_c=corrhat_c,phat0_e1=phat0_e1,phat0_e2=phat0_e2))
}  


