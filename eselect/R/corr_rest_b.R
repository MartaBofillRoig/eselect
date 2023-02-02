##############################################################
#' corr_rest
#' @description computes the correlation restrictions
#' phat0_e1,phat0_e2: estimated probabilities
#' p0_e1,p0_e2,p1_e1,p1_e2: assumed probabilities (design)
#' OR1,OR2: expected effect sizes
#' @return Correlation bounds based on estimated and assumed probabilities
#' @keywords internal
#' @export

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
