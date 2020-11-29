#' eselectsim_ub
#'
#' @description
#'
#' @param ss_arm numeric parameter, sample size per arm
#' @param p0_e1 numeric parameter, probability of occurrence E1 in the control group
#' @param p0_e2 numeric parameter, probability of occurrence E2 in the control group
#' @param OR1 numeric parameter, Odds ratio for the endpoint 1
#' @param OR2 numeric parameter, Odds ratio for the endpoint 2
#' @param p_init numeric parameter, percentage of sample size used in the interim
#' @param criteria decision criteria to choose between the composite endpoint or the endpoint 1 as primary endpoint ("SS": Ratio sample sizes, "ARE": Asymptotic Relative Efficiency).
#'
#' @export
#'
#' @return
#'
#' @details
#'
#' @references
#'

##############################################################
##############################################################
# eselection_b:
# simulation
# estimation correlation blinded without previous info corr/CE
# computation sample size SS
# computation statistic according to the decision (based on SS)

eselectsim_ub <- function(ss_arm,p0_e1,OR1,p0_e2,OR2,p0_ce,p_init=1,criteria="SS",H0_e1=FALSE,H0_e2=FALSE){

  n_init=ss_arm
  ss_arm=round(n_init*p_init)
  total_ss = ss_arm*2

  if(H0_e1==FALSE){
    p1_e1 = (OR1*p0_e1/(1-p0_e1))/(1+(OR1*p0_e1/(1-p0_e1)))
  }else{
    p1_e1=p0_e1
  }

  if(H0_e2==FALSE){
    p1_e2 = (OR2*p0_e2/(1-p0_e2))/(1+(OR2*p0_e2/(1-p0_e2)))

  }else{
    p1_e2=p0_e2
  }

  rho_ds = -(p0_ce-(1-(1-p0_e1)*(1-p0_e2)))/(sqrt(p0_e1*p0_e2*(1-p0_e1)*(1-p0_e2)))
  p1_ce = prob_cbe(p_e1=p1_e1, p_e2=p1_e2, rho=rho_ds)

  # control group
  sm0 = f_sim(samplesize=ss_arm,p_e1=p0_e1,p_e2=p0_e2,p_ce=p0_ce)

  # intervention group
  sm1 = f_sim(samplesize=ss_arm,p_e1=p1_e1,p_e2=p1_e2,p_ce=p1_ce)

  # pooled sample
  sm = sm0 + sm1

  selection = eselect_ub(db0=sm0,db1=sm1,p0_e1=p0_e1,OR1=OR1,p0_e2=p0_e2,OR2=OR2,criteria=criteria)

  if(selection$Decision == 1){

    if(total_ss<selection$SampleSize){
      ss_arm_old = ss_arm
      ss_arm = selection$SampleSize/2

      # if(ss_arm - n_init>0){
        # control group
        sm0_add = f_sim(samplesize=round(ss_arm - ss_arm_old),p_e1=p0_e1,p_e2=p0_e2,p_ce=p0_ce)
        sm0 = sm0 + sm0_add
        # intervention group
        sm1_add = f_sim(samplesize=round(ss_arm - ss_arm_old),p_e1=p1_e1,p_e2=p1_e2,p_ce=p1_ce)
        sm1 = sm1 + sm1_add
      # }
    }

    phat_group1 = 1-(sm1[4])/ss_arm
    phat_group0 = 1-(sm0[4])/ss_arm

    # test odds ratio CE
    TestOR_unpooled = test_f(OR=(phat_group1/(1-phat_group1))/(phat_group0/(1-phat_group0)),p0=phat_group0,n=ss_arm)

  }else{

    if(total_ss<selection$SampleSize){
      ss_arm_old = ss_arm
      ss_arm = selection$SampleSize/2
      # control group
      sm0_add = f_sim(samplesize=round(ss_arm - ss_arm_old),p_e1=p0_e1,p_e2=p0_e2,p_ce=p0_ce)
      sm0 = sm0 + sm0_add
      # intervention group
      sm1_add = f_sim(samplesize=round(ss_arm - ss_arm_old),p_e1=p1_e1,p_e2=p1_e2,p_ce=p1_ce)
      sm1 = sm1 + sm1_add
    }

    phat_group1 = (sm1[1]+sm1[2])/ss_arm
    phat_group0 = (sm0[1]+sm0[2])/ss_arm

    # test odds ratio RE
    TestOR_unpooled = test_f(OR=(phat_group1/(1-phat_group1))/(phat_group0/(1-phat_group0)),p0=phat_group0,n=ss_arm)
  }

  output = list(Test=TestOR_unpooled, Decision=selection$Decision)

  return(output)
}
