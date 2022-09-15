#' Blinded estimation of the correlation
#'
#' @description  This function estimates the correlation between two binary endpoints based on blinded data.
#' @param samplesize numeric parameter, sample size per arm
#' @param p0_e1 numeric parameter, probability of occurrence E1 in the control group
#' @param p0_e2 numeric parameter, probability of occurrence E2 in the control group
#' @param p0_ce numeric parameter, probability of occurrence composite endpoint in the control group
#' @param p1_e1 numeric parameter, probability of occurrence E1 in the treatment group
#' @param p1_e2 numeric parameter, probability of occurrence E2 in the treatment group
#' @param p1_ce numeric parameter, probability of occurrence composite endpoint in the treatment group
#' @param OR1 numeric parameter, Odds ratio for the endpoint 1 (design)
#' @param OR2 numeric parameter, Odds ratio for the endpoint 2 (design)
#' @return This function returns the estimated correlation and the truncated correlation within the possible margins.
#' @keywords internal
#' @export

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

