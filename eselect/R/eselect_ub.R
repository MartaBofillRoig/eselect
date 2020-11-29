#' eselect_ub: Endpoint selection and sample size reassessment for composite endpoints based on unblinded data
#'
#' @description
#'
#' @param db0 matrix
#' @param db1 matrix
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
# eselect_ub_int:
# simulation
# estimation correlation unblinded without previous info corr/CE
# computation ratio sample size
# computation statistic according to the decision

eselect_ub <- function(db0,db1,p0_e1,OR1,p0_e2,OR2,criteria="SS",H0=FALSE){

  samplesize = sum(db0)
  total_ss = 2*samplesize

  p1_e1 = (OR1*p0_e1/(1-p0_e1))/(1+(OR1*p0_e1/(1-p0_e1)))
  p1_e2 = (OR2*p0_e2/(1-p0_e2))/(1+(OR2*p0_e2/(1-p0_e2)))

  # estimated probabilities
  phat0_e1 = (db0[1]+db0[2])/samplesize
  phat0_e2 = (db0[1]+db0[3])/samplesize
  phat0_ce = 1-(db0[4])/samplesize

  phat1_e1 = (db1[1]+db1[2])/samplesize
  phat1_e2 = (db1[1]+db1[3])/samplesize
  phat1_ce = 1-(db1[4])/samplesize

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

  if(criteria=="SS"){

    ratio_ss = samplesize_e1/samplesize_estar

    if(ratio_ss>= 1){

      total_ss=samplesize_estar
      Decision = 1

    }else{

      total_ss=samplesize_e1
      Decision = 0

    }
  }
  if(criteria=="ARE"){

    ARE_up = ARE_cbe(p0_e1=phat0_e1, p0_e2=phat0_e2, eff_e1=OR1, eff_e2=OR2, rho=corrhat_c)

    if(ARE_up>= 1){

      total_ss=samplesize_estar
      Decision = 1

    }else{

      total_ss=samplesize_e1
      Decision = 0

    }
  }

  return(list(SampleSize = total_ss, Decision = Decision))
}
