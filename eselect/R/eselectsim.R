#' Simulation trials with endpoint selection and sample size reassessment for composite endpoints based on blinded data
#'
#' @description  This function simulates trials with endpoint selection and sample size reassessment  for composite binary endpoints based on blinded data. The composite endpoint is assumed to be a binary endpoint formed by a combination of two events (E1 and E2). We assume that the endpoint 1 is more relevant for the clinical question than endpoint 2. This function simulates a trial based on the design parameters and use the algorithm implemented in eselect() to select the primary endpoint and recalculate the sample size accordingly.
#'
#' @param ss_arm numeric parameter, sample size per arm
#' @param p0_e1 numeric parameter, probability of occurrence E1 in the control group
#' @param p0_e2 numeric parameter, probability of occurrence E2 in the control group
#' @param p0_ce numeric parameter, probability of occurrence composite endpoint in the control group
#' @param OR1 numeric parameter, Odds ratio for the endpoint 1
#' @param OR2 numeric parameter, Odds ratio for the endpoint 2
#' @param p_init numeric parameter, percentage of sample size used in the interim
#' @param criteria decision criteria to choose between the composite endpoint or the endpoint 1 as primary endpoint ("SS": Ratio sample sizes, "ARE": Asymptotic Relative Efficiency).
#' @param H0_e1 Simulate under true null hypothesis for the endpoint E1 (TRUE/FALSE).
#' @param H0_e2 Simulate under true null hypothesis for the endpoint E2 (TRUE/FALSE).
#' @param SS_r Sample size reassessment (TRUE/FALSE). If TRUE, in those cases where the sample size is less than the needed for achieving the pre-specified power, additional subjects are added after recalculating the sample size. If FALSE, no more subjects are added in the study.
#' @param alpha Type I error.
#' @param beta Type II error.
#'
#' @export
#' @import CompAREdesign
#' @references Bofill Roig, M., Gómez Melis, G., Posch, M., & Koenig, F. (2022). Adaptive clinical trial designs with blinded selection of binary composite endpoints and sample size reassessment. Biostatistics (in press). arXiv e-prints, arXiv-2206 (https://doi.org/10.48550/arXiv.2206.09639).
#' @return This function returns the decision (Decision = 1, meaning the chosen endpoint is the composite endpoint; and Decision = 0, meaning the chosen endpoint is the relevant endpoint) and the statistic to test the primary hypothesis according to the decision.
#'
#'
#'

##############################################################
##############################################################
# eselection_b:
# simulation
# estimation correlation blinded without previous info corr/CE
# computation sample size SS
# computation statistic according to the decision (based on SS)

eselectsim <- function(ss_arm,p0_e1,OR1,p0_e2,OR2,p0_ce,p_init=1,criteria="SS",H0_e1=FALSE,H0_e2=FALSE,SS_r=TRUE,alpha=0.05,beta=0.2){

  requireNamespace("CompAREdesign")

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

  if(rho_ds<lower_corr(p_e1=p1_e1, p_e2=p1_e2)){
    rho_ds=lower_corr(p_e1=p1_e1, p_e2=p1_e2)
  }
  if(rho_ds>upper_corr(p_e1=p1_e1, p_e2=p1_e2)){
    rho_ds=upper_corr(p_e1=p1_e1, p_e2=p1_e2)
  }

  p1_ce = prob_cbe(p_e1=p1_e1, p_e2=p1_e2, rho=rho_ds)

  # control group
  sm0 = f_sim(samplesize=ss_arm,p_e1=p0_e1,p_e2=p0_e2,p_ce=p0_ce)

  # intervention group
  sm1 = f_sim(samplesize=ss_arm,p_e1=p1_e1,p_e2=p1_e2,p_ce=p1_ce)

  # pooled sample
  sm = sm0 + sm1

  selection = eselect(db=sm,p0_e1=p0_e1,OR1=OR1,p0_e2=p0_e2,OR2=OR2,criteria=criteria,alpha=alpha,beta=beta)

  if(selection$Decision == 1){

    if(total_ss<selection$SampleSize & SS_r==TRUE){
      ss_arm_old = ss_arm
      ss_arm = selection$SampleSize/2

      # control group
      sm0_add = f_sim(samplesize=round(ss_arm - ss_arm_old),p_e1=p0_e1,p_e2=p0_e2,p_ce=p0_ce)
      sm0 = sm0 + sm0_add
      # intervention group
      sm1_add = f_sim(samplesize=round(ss_arm - ss_arm_old),p_e1=p1_e1,p_e2=p1_e2,p_ce=p1_ce)
      sm1 = sm1 + sm1_add
    }
    if(total_ss<2*n_init & SS_r==FALSE){
      ss_arm_old = ss_arm
      ss_arm = n_init

      # control group
      sm0_add = f_sim(samplesize=round(ss_arm - ss_arm_old),p_e1=p0_e1,p_e2=p0_e2,p_ce=p0_ce)
      sm0 = sm0 + sm0_add
      # intervention group
      sm1_add = f_sim(samplesize=round(ss_arm - ss_arm_old),p_e1=p1_e1,p_e2=p1_e2,p_ce=p1_ce)
      sm1 = sm1 + sm1_add
    }

    phat_group1 = 1-(sm1[4])/ss_arm
    phat_group0 = 1-(sm0[4])/ss_arm

    # test odds ratio CE
    TestOR_unpooled = test_f(OR=(phat_group1/(1-phat_group1))/(phat_group0/(1-phat_group0)),p0=phat_group0,n=ss_arm)

  }else{

    if(total_ss<selection$SampleSize & SS_r==TRUE){
      ss_arm_old = ss_arm
      ss_arm = selection$SampleSize/2
      # control group
      sm0_add = f_sim(samplesize=round(ss_arm - ss_arm_old),p_e1=p0_e1,p_e2=p0_e2,p_ce=p0_ce)
      sm0 = sm0 + sm0_add
      # intervention group
      sm1_add = f_sim(samplesize=round(ss_arm - ss_arm_old),p_e1=p1_e1,p_e2=p1_e2,p_ce=p1_ce)
      sm1 = sm1 + sm1_add
    }
    if(total_ss<2*n_init & SS_r==FALSE){
      ss_arm_old = ss_arm
      ss_arm = n_init

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
