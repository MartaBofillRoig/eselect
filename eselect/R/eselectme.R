#' eselectme: Endpoint selection and sample size reassessment for composite endpoints based on blinded data
#'
#' @description
#'
#' @param db matrix 2x2 table (pooled sample)
#' @param p0_e1 numeric parameter, probability of occurrence E1 in the control group
#' @param p0_e2 numeric parameter, probability of occurrence E2 in the control group
#' @param OR1 numeric parameter, Odds ratio for the endpoint 1
#' @param OR2 numeric parameter, Odds ratio for the endpoint 2
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

eselectme <- function(db,p0_e1,OR1,p0_e2,OR2){

  total_ss = sum(db)
  ss_arm = total_ss/2

  p1_e1 = (OR1*p0_e1/(1-p0_e1))/(1+(OR1*p0_e1/(1-p0_e1)))
  p1_e2 = (OR2*p0_e2/(1-p0_e2))/(1+(OR2*p0_e2/(1-p0_e2)))

  # estimated probabilities
  phat_e1 = (db[1]+db[2])/total_ss
  phat_e2 = (db[1]+db[3])/total_ss
  phat_ce = 1-(db[4])/total_ss

  #
  phat0_e1 = fun_p0(p=phat_e1,l=OR1)
  phat0_e2 = fun_p0(p=phat_e2,l=OR2)

  #
  phat1_e1 = (OR1*phat0_e1/(1-phat0_e1))/(1+(OR1*phat0_e1/(1-phat0_e1)))
  phat1_e2 = (OR2*phat0_e2/(1-phat0_e2))/(1+(OR2*phat0_e2/(1-phat0_e2)))

  # estimated correlation
  corrhat =  (phat_ce - (ss_arm/total_ss)*(1 - (1-phat0_e1)*(1-phat0_e2)) - (ss_arm/total_ss)*(1-(1-phat1_e1)*(1-phat1_e2)))/(-(ss_arm/total_ss)*sqrt(phat0_e1*phat0_e2*(1-phat0_e1)*(1-phat0_e2))-(ss_arm/total_ss)*sqrt(phat1_e1*phat1_e2*(1-phat1_e1)*(1-phat1_e2)))

  # correlation restrictions
  rest = corr_rest_b(phat0_e1,phat0_e2,p0_e1,p0_e2,p1_e1,p1_e2,OR1,OR2)
  low = rest[1]
  upp = rest[2]

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

  # criteria=="SS"

  ratio_ss = samplesize_e1/samplesize_estar

  if(ratio_ss>= 1){

    total_ss=samplesize_estar
    Decision = 1

    # Shoud we use CE and RE as multiple primary endpoints?
    # Comparison (CE,RE) vs RE - Simulation-based

    alpha= 0.05
    z.alpha <- qnorm(1-alpha,0,1)

    # simulation multiple test (adjusted bonferroni)
    corr_struct=diag(2)
    corr_struct[1,2]=corrhat_c
    corr_struct[2,1]=corrhat_c
    quant <- qmvnorm(1-alpha, corr = corr_struct, tail = "lower.tail")

    # simtest_me(OR1=OR2,p0_e1=phat0_e1,OR2=OR2,p0_e2=phat0_e2,n=ss_arm)< -quant$quantile
    power_multendpoints <- sum(replicate(nsim,(simtest_me(OR1=OR2,p0_e1=p0_e1,OR2=OR2,p0_e2=p0_e2,n=ss_arm)< -quant$quantile)>0))/nsim

    # simulation composite test
    power_CE <- sum(replicate(nsim,f_OR(samplesize=ss_arm,
                                        p0=p0_ce,
                                        p1=prob_cbe(p_e1=p1_e1, p_e2=p1_e2, rho=corr_struct)))< - z.alpha)/nsim

    if(power_CE<power_multendpoints){
      Decision = 2
    }

  }else{

    total_ss=samplesize_e1
    Decision = 0

  }

  return(list(SampleSize = total_ss, Decision = Decision))
}


