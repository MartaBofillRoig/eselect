#' Unblinded estimation of the correlation
#'
#' @description  This function estimates the correlation between two binary endpoints based on unblinded data.
#' @param samplesize numeric parameter, sample size per arm
#' @param p0_e1 numeric parameter, probability of occurrence E1 in the control group
#' @param p0_e2 numeric parameter, probability of occurrence E2 in the control group
#' @param p0_ce numeric parameter, probability of occurrence composite endpoint in the control group
#' @param p1_e1 numeric parameter, probability of occurrence E1 in the treatment group
#' @param p1_e2 numeric parameter, probability of occurrence E2 in the treatment group
#' @param p1_ce numeric parameter, probability of occurrence composite endpoint in the treatment group
#' @param OR1 numeric parameter, Odds ratio for the endpoint 1 (design)
#' @param OR2 numeric parameter, Odds ratio for the endpoint 2 (design)
#' @return This function returns the estimated correlation based on unblinded data and the truncated correlation within the possible margins.
#' @keywords internal
#' @export
#'
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
