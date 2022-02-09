##############################################################
#' simtest_me
#' @description simulates two correlated binary endpoints and computes the statistical test (OR) for two endpoints
#' @keywords internal
#'
#'

simtest_me <- function(OR1,p0_e1,OR2,p0_e2,p0_ce,p1_ce,n,ss_arm){
  
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