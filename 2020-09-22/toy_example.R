

##############################################################  
i=1
p0_e1 = dataset$p0_e1[i]
p1_e1 = dataset$p1_e1[i]
p0_e2 = dataset$p0_e2[i]
p1_e2 = dataset$p1_e2[i]
p0_ce = dataset$p0_ce[i]
p1_ce = dataset$p1_ce[i]
OR1 =dataset$OR1[i]
samplesize=500
# dataset$corr[i]
OR2 =dataset$OR2[i]
OR_ce=dataset$OR_ce[i]

upp =dataset$max_corr[i]
low =dataset$min_corr[i]
dataset[i,]

# control group
s1_group0 = p0_e1 + p0_e2 - p0_ce
s2_group0 = p0_ce-p0_e2  
s3_group0 = p0_ce-p0_e1  
s4_group0 = 1- p0_ce  

sm0 = rmultinom(1,500,c(s1_group0,s2_group0,s3_group0,s4_group0))
corr_est0 =(sm0[1]*sm0[4]-(sm0[2]*sm0[3]))/sqrt((sm0[1]+sm0[3])*(sm0[2]+sm0[4])*(sm0[1]+sm0[2])*(sm0[3]+sm0[4]))
(corr_est0)

# intervention group
s1_group1 = p1_e1+p1_e2-p1_ce
s2_group1 = ifelse(p1_ce-p1_e2>0, p1_ce-p1_e2, 0)  
s3_group1 = ifelse(p1_ce-p1_e1>0, p1_ce-p1_e1, 0)   
s4_group1 = 1- p1_ce  

sm1 = rmultinom(1,500,c(s1_group1,s2_group1,s3_group1,s4_group1))
corr_est1 =(sm1[1]*sm1[4]-(sm1[2]*sm1[3]))/sqrt((sm1[1]+sm1[3])*(sm1[2]+sm1[4])*(sm1[1]+sm1[2])*(sm1[3]+sm1[4]))
(corr_est1)


sm = (sm0+sm1)
corr_est =(sm[1]*sm[4]-(sm[2]*sm[3]))/sqrt((sm[1]+sm[3])*(sm[2]+sm[4])*(sm[1]+sm[2])*(sm[3]+sm[4]))
(corr_est) 

##############################################################  

dataset$OR_ce = (dataset$p1_ce/(1-dataset$p1_ce))/(dataset$p0_ce/(1-dataset$p0_ce))

##############################################################  

##############################################################  

# f_ES <- function(samplesize,p0_e1,p1_e1,p0_e2,p1_e2,p0_ce,p1_ce,upp,low){ 
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
  
  phat0_e1 = fun_p0(p=phat_e1,l=OR1)
  phat0_e2 = fun_p0(p=phat_e2,l=OR2)
  phat0_ce = fun_p0(p=phat_ce,l=OR_ce)
  
  corrhat = ((phat0_e1+phat0_e2-phat0_ce)-phat0_e1*phat0_e2)/sqrt(phat0_e1*(1-phat0_e1)*phat0_e2*(1-phat0_e2))
  
  if(corrhat > upp){
    corrhat_c = upp
  }else if(corrhat < low){
    corrhat_c = low
  }else if(corrhat >= low && corrhat <= upp){
    corrhat_c = corrhat
  }
  
  # decision based on blinded data
  ARE_up = ARE_cbe(p0_e1=phat0_e1, p0_e2=phat0_e2, eff_e1=OR1, eff_e2=OR2, rho=corrhat_c) 
  
  if(ARE_up>= 1){
    phat_group1 = 1-(sm1[4])/ss
    phat_group0 = 1-(sm0[4])/ss
    
    # test odds ratio RE
    TestOR_unpooled = log((phat_group1/(1-phat_group1))/(phat_group0/(1-phat_group0)))*((1/(phat_group0*(1-phat_group0))+ 1/(phat_group1*(1-phat_group1)))/samplesize)^(-1/2)
  }else{
    phat_group1 = (sm1[1]+sm1[2])/ss
    phat_group1 = (sm0[1]+sm0[2])/ss
    
    # test odds ratio RE
    TestOR_unpooled = log((phat_group1/(1-phat_group1))/(phat_group0/(1-phat_group0)))*((1/(phat_group0*(1-phat_group0))+ 1/(phat_group1*(1-phat_group1)))/samplesize)^(-1/2)
  }
  
  # return(TestOR_unpooled)
# }