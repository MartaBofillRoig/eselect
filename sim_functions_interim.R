

##############################################################  
##############################################################  
# eselection_b:
# simulation
# estimation correlation blinded without previous info corr/CE
# computation sample size SS
# computation statistic according to the decision (based on SS)

eselection_bSS <- function(samplesize,p0_e1,p1_e1,OR1,p0_e2,p1_e2,OR2,p0_ce,p1_ce){  
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
  
  # decision based on blinded data
  samplesize_e1 = samplesize_OR(p0=phat0_e1, OR=OR1, alpha=0.05, beta=0.2)
  samplesize_estar = samplesize_cbe(p0_e1=phat0_e1,p0_e2=phat0_e2,
                                    eff_e1 = OR1,
                                    effm_e1 = "or",
                                    eff_e2 = OR2,
                                    effm_e2 = "or",
                                    effm_ce = "or",
                                    rho=corrhat_c,
                                    alpha = 0.05,
                                    beta = 0.2)
  
  ARE_up = samplesize_e1/samplesize_estar
  
  if(ARE_up>= 1){
    phat_group1 = 1-(sm1[4])/samplesize
    phat_group0 = 1-(sm0[4])/samplesize
    
    # test odds ratio CE
    TestOR_unpooled = log((phat_group1/(1-phat_group1))/(phat_group0/(1-phat_group0)))*((1/(phat_group0*(1-phat_group0))+ 1/(phat_group1*(1-phat_group1)))/samplesize)^(-1/2)
    
    Decision = 1
    
  }else{
    phat_group1 = (sm1[1]+sm1[2])/samplesize
    phat_group0 = (sm0[1]+sm0[2])/samplesize
    
    # test odds ratio RE
    TestOR_unpooled = log((phat_group1/(1-phat_group1))/(phat_group0/(1-phat_group0)))*((1/(phat_group0*(1-phat_group0))+ 1/(phat_group1*(1-phat_group1)))/samplesize)^(-1/2)
    
    Decision = 0
  }
  
  return(c(TestOR_unpooled,Decision))
}


##############################################################  
# eselection_ub:
# simulation
# estimation correlation unblinded without previous info corr/CE
# computation ratio sample size
# computation statistic according to the decision

eselection_ub_int <- function(samplesize,p0_e1,p1_e1,OR1,p0_e2,p1_e2,OR2,p0_ce,p1_ce,p_init=1){ 
  
  n_init=samplesize
  samplesize=round(n_init*p_init)
  
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
  
  rest = corr_rest_ub(phat0_e1=phat0_e1,phat0_e2=phat0_e2,
                      phat1_e1=phat1_e1,phat1_e2=phat1_e2,
                      p0_e1=p0_e1,p0_e2=p0_e2,
                      p1_e1=p1_e1,p1_e2=p1_e2,
                      OR1=OR1,OR2=OR2)
  
  upp = rest[2]
  low = rest[1] 
  
  corrhat_c = corrhat
  if(corrhat > upp){
    corrhat_c = upp
  }
  if(corrhat < low){
    corrhat_c = low
  } 
  
  samplesize_e1 = samplesize_OR(p0=phat0_e1, OR=OR1, alpha=0.05, beta=0.2)
  samplesize_estar = samplesize_cbe(p0_e1=phat0_e1,p0_e2=phat0_e2,
                                    eff_e1 = OR1,
                                    effm_e1 = "or",
                                    eff_e2 = OR2,
                                    effm_e2 = "or",
                                    effm_ce = "or",
                                    rho=corrhat_c,
                                    alpha = 0.05,
                                    beta = 0.2)
  
  ratio_ss = samplesize_e1/samplesize_estar
  
  if(ratio_ss>= 1){ 
    
    if(total_ss<samplesize_estar){
      # control group
      sm0_add = f_sim(samplesize=round(samplesize_estar/2 - samplesize),p_e1=p0_e1,p_e2=p0_e2,p_ce=p0_ce)
      sm0 = sm0 + sm0_add
      # intervention group
      sm1_add = f_sim(samplesize=round(samplesize_estar/2 - samplesize),p_e1=p1_e1,p_e2=p1_e2,p_ce=p1_ce)
      sm1 = sm1 + sm1_add  
    }
    
    phat_group1 = 1-(sm1[4])/samplesize
    phat_group0 = 1-(sm0[4])/samplesize
    
    # test odds ratio CE
    TestOR_unpooled = test_f(OR=(phat_group1/(1-phat_group1))/(phat_group0/(1-phat_group0)),p0=phat_group0,n=samplesize) 
    
    Decision = 1
    
  }else{
    
    if(total_ss<samplesize_e1){
      # control group
      sm0_add = f_sim(samplesize=round(samplesize_e1/2 - samplesize),p_e1=p0_e1,p_e2=p0_e2,p_ce=p0_ce)
      sm0 = sm0 + sm0_add
      # intervention group
      sm1_add = f_sim(samplesize=round(samplesize_e1/2 - samplesize),p_e1=p1_e1,p_e2=p1_e2,p_ce=p1_ce)
      sm1 = sm1 + sm1_add  
    }
    
    phat_group1 = (sm1[1]+sm1[2])/samplesize
    phat_group0 = (sm0[1]+sm0[2])/samplesize
    
    # test odds ratio RE
    TestOR_unpooled = test_f(OR=(phat_group1/(1-phat_group1))/(phat_group0/(1-phat_group0)),p0=phat_group0,n=samplesize) 
    
    Decision = 0
  }
  
  return(c(TestOR_unpooled,Decision))
}
