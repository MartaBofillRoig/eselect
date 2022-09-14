#' Endpoint selection and sample size reassessment for composite endpoints based on blinded data
#'
#' @description Endpoint selection and sample size reassessment for composite endpoints based on blinded data. The composite endpoint is assumed to be a binary endpoint formed by a combination of two events (E1 and E2). We assume that the endpoint 1 is more relevant for the clinical question than endpoint 2. This function selects between the composite endpoint or the relevant endpoint as the primary endpoint of the study and recalculate the sample size accordingly. The decision criteria to decide between the composite endpoint or the relevant endpoint might be the ratio of the corresponding sample sizes ("SS") or the Asymptotic Relative Efficiency ("ARE").
#' The algorithm of the function is the following: First, the probabilities of the composite components in the control group and the correlation between them are estimated based on blinded data. Second, using the estimated probabilities and the estimated correlation, the decision criteria is computed and the primary endpoint is selected.  Finally, the sample size is recalculated according to the decision.
#'
#' @param db matrix 2x2 table (pooled sample)
#' @param p0_e1 numeric parameter, probability of occurrence E1 in the control group
#' @param p0_e2 numeric parameter, probability of occurrence E2 in the control group
#' @param OR1 numeric parameter, Odds ratio for the endpoint 1
#' @param OR2 numeric parameter, Odds ratio for the endpoint 2
#' @param criteria decision criteria to choose between the composite endpoint or the endpoint 1 as primary endpoint ("SS": Ratio sample sizes, "ARE": Asymptotic Relative Efficiency).
#' @param alpha Type I error.
#' @param beta Type II error.
#'
#' @export
#' @import CompAREdesign
#' @references Bofill Roig, M., Gómez Melis, G., Posch, M., & Koenig, F. (2022). Adaptive clinical trial designs with blinded selection of binary composite endpoints and sample size reassessment. Biostatistics (in press). arXiv e-prints, arXiv-2206 (https://doi.org/10.48550/arXiv.2206.09639).
#' Bofill Roig, M., & Gómez Melis, G. "Selection of composite binary endpoints in clinical trials." Biometrical Journal 60.2 (2018): 246-261.
#' @return This function returns the decision (Decision = 1, meaning the chosen endpoint is the composite endpoint; and Decision = 0, meaning the chosen endpoint is the relevant endpoint) and the sample size according to the decision.
#'
#' @examples
#' # Based on Bofill Roig, M., et al.
#' # (See supplementary material in https://doi.org/10.48550/arXiv.2206.09639)
#' p0_e1 = 0.173
#' p0_e2 = 0.055
#' p1_e1 = 0.121;
#' p1_e2 = 0.057;
#' n1 = 569
#' n0 = 576
#' n = n0+n1
#' p1 = (p0_e1*n0 + p1_e1*n1)/n
#' p2 = (p0_e2*n0 + p1_e2*n1)/n
#' p_ce = (0.203*n0 + 0.146*n1)/n
#' OR1 = 0.7
#' OR2 = 0.9
#' x11 = round((p1+p2-p_ce)*n)
#' x12 = round((p1)*n-x11)
#' x21 = round((p2)*n- x11)
#' x22 = round((1-p_ce)*n)
#' data = matrix(c(x11,x12,x21,x22), nrow = 2 , ncol = 2, byrow = FALSE)
#' eselect(db=data,p0_e1=0.18,OR1=0.70,p0_e2=0.05,OR2=0.9,criteria="SS",alpha=0.05,beta=0.2)
#'
#'

##############################################################
##############################################################
# eselection_b:

eselect <- function(db,p0_e1,OR1,p0_e2,OR2,criteria="SS",alpha=0.05,beta=0.2){

  requireNamespace("CompAREdesign")

  # if(is.table(db)==F || (dim(db)==c(3,3))==c(F,F) || is.numeric(db)==F){
  #   stop("The table must be a matrix 2x2 table")
  # }else
  if(p0_e1 < 0 || p0_e1 > 1){
    stop("The probability of observing the event E1 (p_e1) must be number between 0 and 1")
  }else if(p0_e2 < 0 || p0_e2 > 1){
    stop("The probability of observing the event E2 (p_e2) must be number between 0 and 1")
  }else if(OR1 < 0 || OR1 >= 1){
    stop("The odds ratio (OR1) must be number between 0 and 1")
  }else if(OR2 < 0 || OR2 >= 1){
    stop("The odds ratio (OR2) must be number between 0 and 1")
  }else if(criteria != "SS" && criteria != "ARE"){
    stop("You have to choose between sample size (SS) or Asymptotic Relative Efficiency (ARE)")
  }else if( 0 > alpha || alpha > 1){
    stop("Alpha value must be number between 0 and 1")
  }else if( 0 > beta || beta > 1){
    stop("Beta value must be number between 0 and 1")
  }

  #

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

  samplesize_e1 = samplesize_OR(p0=phat0_e1, OR=OR1, alpha=alpha, beta=beta)
  samplesize_estar = samplesize_cbe(p0_e1=phat0_e1,p0_e2=phat0_e2,
                                    eff_e1 = OR1,
                                    effm_e1 = "or",
                                    eff_e2 = OR2,
                                    effm_e2 = "or",
                                    effm_ce = "or",
                                    rho=corrhat_c,
                                    alpha = alpha,
                                    beta = beta)

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


