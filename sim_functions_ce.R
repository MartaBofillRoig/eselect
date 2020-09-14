

##############################################################  

OR_function<- function(p0, p1){
  
  OR<- (p1/(1-p1))/(p0/(1-p0))
  return(OR)
}

##############################################################  

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

# p=p1_e1+p0_e1
# l=OR1
# fun_p0(p=p,l=l)

##############################################################  

f_ES <- function(samplesize,p0_e1,p1_e1,OR1,p0_e2,p1_e2,OR2,p0_ce,p1_ce,OR_ce,upp,low){ 
  ss=samplesize/2
  
  # control group
  s1_group0 = p0_e1 + p0_e2 - p0_ce
  s2_group0 = p0_ce-p0_e2  
  s3_group0 = p0_ce-p0_e1  
  s4_group0 = 1- p0_ce  
  
  sm0 = rmultinom(1,ss,c(s1_group0,s2_group0,s3_group0,s4_group0)) 
  
  # intervention group
  s1_group1 = p1_e1+p1_e2-p1_ce
  s2_group1 = ifelse(p1_ce-p1_e2>0, p1_ce-p1_e2, 0)  
  s3_group1 = ifelse(p1_ce-p1_e1>0, p1_ce-p1_e1, 0)   
  s4_group1 = 1- p1_ce  
  
  sm1 = rmultinom(1,ss,c(s1_group1,s2_group1,s3_group1,s4_group1))
  
  # pooled sample
  sm = sm0 + sm1
  
  phat_e1 = (sm[1]+sm[2])/samplesize
  phat_e2 = (sm[1]+sm[3])/samplesize
  phat_ce = 1-(sm[4])/samplesize
  
  # OR1 = OR_function(p0_e1,p1_e1)
  # OR2 = OR_function(p0_e2,p1_e2)
  # OR_ce = OR_function(p0_ce,p1_ce)
    
  phat0_e1 = fun_p0(p=phat_e1,l=OR1)
  phat0_e2 = fun_p0(p=phat_e2,l=OR2)
  phat0_ce = fun_p0(p=phat_ce,l=OR_ce)
  
  corrhat = ((phat0_e1+phat0_e2-phat0_ce)-phat0_e1*phat0_e2)/sqrt(phat0_e1*(1-phat0_e1)*phat0_e2*(1-phat0_e2))
  
  update_uppcorr0=upper_corr(phat0_e1,phat0_e2)
  update_uppcorr1=upper_corr((OR1*phat0_e1/(1-phat0_e1))/(1+(OR1*phat0_e1/(1-phat0_e1))),(OR2*phat0_e2/(1-phat0_e2))/(1+(OR2*phat0_e2/(1-phat0_e2))))
    
  update_lowcorr0= lower_corr(phat0_e1,phat0_e2)  
  update_lowcorr1=lower_corr((OR1*phat0_e1/(1-phat0_e1))/(1+(OR1*phat0_e1/(1-phat0_e1))),(OR2*phat0_e2/(1-phat0_e2))/(1+(OR2*phat0_e2/(1-phat0_e2))))
  
  corrhat_c = corrhat
  if(corrhat > upp || corrhat > update_uppcorr0 || corrhat > update_uppcorr1){
    corrhat_c = min(upp,upper_corr(phat0_e1,phat0_e2))
  }
  if(corrhat < low || corrhat < update_lowcorr0 || corrhat < update_lowcorr1){
    corrhat_c = max(low,lower_corr(phat0_e1,phat0_e2),lower_corr((OR1*phat0_e1/(1-phat0_e1))/(1+(OR1*phat0_e1/(1-phat0_e1))),(OR2*phat0_e2/(1-phat0_e2))/(1+(OR2*phat0_e2/(1-phat0_e2)))))
  } 
  
  # decision based on blinded data
  ARE_up = ARE_cbe(p0_e1=phat0_e1, p0_e2=phat0_e2, eff_e1=OR1, eff_e2=OR2, rho=corrhat_c) 
  # p0_e1=phat0_e1; p0_e2=phat0_e2; eff_e1=OR1; eff_e2=OR2; rho=corrhat_c
  
  if(ARE_up>= 1){
    phat_group1 = 1-(sm1[4])/ss
    phat_group0 = 1-(sm0[4])/ss
    
    # test odds ratio RE
    TestOR_unpooled = log((phat_group1/(1-phat_group1))/(phat_group0/(1-phat_group0)))*((1/(phat_group0*(1-phat_group0))+ 1/(phat_group1*(1-phat_group1)))/ss)^(-1/2)
  }else{
    phat_group1 = (sm1[1]+sm1[2])/ss
    phat_group0 = (sm0[1]+sm0[2])/ss
    
    # test odds ratio RE
    TestOR_unpooled = log((phat_group1/(1-phat_group1))/(phat_group0/(1-phat_group0)))*((1/(phat_group0*(1-phat_group0))+ 1/(phat_group1*(1-phat_group1)))/ss)^(-1/2)
  }
  
  return(TestOR_unpooled)
}