##############################################################
#' Simulation 2x2 table binary endpoints
#' @description simulation binary data two outcomes
#' @param samplesize sample size simulated sample
#' @param p_e1 numeric parameter, probability of occurrence E1
#' @param p_e2 numeric parameter, probability of occurrence E2
#' @param p_ce numeric parameter, probability of occurrence composite endpoint
#' @keywords internal
#' @return Simlated data
#' @export
#' @return the function returns 2x2 table
#' s1+s2: num X1
#' s1+s3: num X2

f_sim <- function(samplesize,p_e1,p_e2,p_ce){

  # 2x2 table
  s1_group = p_e1 + p_e2 - p_ce
  s2_group = ifelse(p_ce-p_e2>0, p_ce-p_e2, 0)#p_ce-p_e2
  s3_group = ifelse(p_ce-p_e1>0, p_ce-p_e1, 0)#p_ce-p_e1
  s4_group = 1- p_ce

  data = rmultinom(n=1,size=round(samplesize),prob=c(s1_group,s2_group,s3_group,s4_group))

  return(data)
}
