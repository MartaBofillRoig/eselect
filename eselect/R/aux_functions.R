
##############################################################
#' samplesize_OR
#' @description Sample size calculations in terms of odds ratio
#' @keywords internal
#' @export
#'

samplesize_OR <- function(p0, OR, alpha=0.05, beta=0.2, Unpooled="Unpooled Variance"){
  z.alpha <- qnorm(1-alpha,0,1)
  z.beta <-  qnorm(1-beta,0,1)

  p1 = (OR*p0/(1-p0))/(1+(OR*p0/(1-p0)))

  if(Unpooled=="Unpooled Variance"){
    # sample size per group
    n1 = ((z.alpha+z.beta)/(log(OR)))^2*( 1/(p0*(1-p0)) + 1/(p1*(1-p1)) )
  }else{
    p = (p1 + p0)/2
    # sample size per group
    n1 = ((z.alpha* sqrt(2/(p*(1-p))) + z.beta* sqrt(1/(p0*(1-p0)) + 1/(p1*(1-p1))))/(log(OR)))^2
  }

  n = 2*n1

  return(n)
}


##############################################################
#' f_OR
#' @description simulates binary outcomes and computes the statistic
#' @keywords internal
#' @export
#' @import stats
#'

f_OR <- function(samplesize,p0,p1){

  # estimation of the probability of observing the composite event in control group
  R = runif(samplesize)
  phat_group0 = 1 - sum(R>=p0 & R<1)/samplesize

  # estimation of the probability of observing the composite event in test group
  R = runif(samplesize)
  phat_group1 = 1 - sum(R>=p1 & R<1)/samplesize

  # test odds ratio
  TestOR_unpooled = test_f(OR=OR_function(p0=phat_group0,p1=phat_group1),p0=phat_group0,n=samplesize)

  return(TestOR_unpooled)
}

##############################################################
#' test_f
#' @description computes the statistical test (OR)
#' @keywords internal
#' @export
#'

test_f <- function(OR,p0,n){

  # OR: odds ratio
  # p0: probability under control group
  # n: sample size per group

  p1 <- (OR*p0/(1-p0))/(1+(OR*p0/(1-p0)))
  den <- 1/(n*(p0*(1-p0))) + 1/(n*(p1*(1-p1)))

  test <- (log(OR))/sqrt(den)
  return(test)

}

##############################################################
#' test_me
#' @description computes the statistical tests (OR) for two endpoints
#' @keywords internal
#' @export
#'

test_me <- function(OR1,p0_e1,OR2,p0_e2,n){

  # OR: (estimated) odds ratio
  # p0: (estimated) probability under control group
  # n: sample size per group

  p1 <- (OR1*p0_e1/(1-p0_e1))/(1+(OR1*p0_e1/(1-p0_e1)))
  den <- 1/(n*(p0_e1*(1-p0_e1))) + 1/(n*(p1*(1-p1)))
  test1 <- (log(OR1))/sqrt(den)

  p1 <- (OR2*p0_e2/(1-p0_e2))/(1+(OR2*p0_e2/(1-p0_e2)))
  den <- 1/(n*(p0_e2*(1-p0_e2))) + 1/(n*(p1*(1-p1)))
  test2 <- (log(OR2))/sqrt(den)

  return(list(Test1=test1,Test2=test2))

}

##############################################################
#' OR_function
#' @description computes the odds ratio
#' @keywords internal
#' @export
#'

OR_function<- function(p0, p1){

  OR<- (p1/(1-p1))/(p0/(1-p0))
  return(OR)
}

##############################################################
#' fun_p0
#' @description computes the probability under the control group
#' based on the pooled probability and the odds ratio
#' p: pooled probability
#' l: odds ratio
#' @keywords internal
#' @export


fun_p0 <- function(p,l){

  p = p*2

  sol1= -(1 + l + p - l*p + sqrt(4*(-1 + l)*p + (1 + l + p - l*p)^2))/(2*(-1 + l))
  sol2= -(1 + l + p - l*p - sqrt(4*(-1 + l)*p + (1 + l + p - l*p)^2))/(2*(-1 + l))

  if(sol2>0 & sol2<1){
    p0=sol2
  }else{
    p0=sol1
  }

  return(p0)
}
# EXAMPLE
# p=p1_e1+p0_e1
# l=OR1
# fun_p0(p=p,l=l)

##############################################################
#' corr_rest_ub: computes the correlation restrictions unblinded case
#' @description computes the correlation restrictions with unblinded data
#' phat0_e1,phat0_e2,phat1_e1,phat1_e2: estimated probabilities
#' p0_e1,p0_e2,p1_e1,p1_e2: assumed probabilities (design)
#' OR1,OR2: expected effect sizes
#' @keywords internal
#' @export

corr_rest_ub <- function(phat0_e1,phat0_e2,phat1_e1,phat1_e2,p0_e1,p0_e2,p1_e1,p1_e2,OR1,OR2){
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

  return(c(low,upp))
}

##############################################################
#' corr_rest
#' @description computes the correlation restrictions
#' phat0_e1,phat0_e2: estimated probabilities
#' p0_e1,p0_e2,p1_e1,p1_e2: assumed probabilities (design)
#' OR1,OR2: expected effect sizes
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


##############################################################
#' Simulation 2x2 table binary endpoints
#' @description simulation binary data two outcomes
#' @param samplesize sample size simulated sample
#' @param p_e1 numeric parameter, probability of occurrence E1
#' @param p_e2 numeric parameter, probability of occurrence E2
#' @param p_ce numeric parameter, probability of occurrence composite endpoint
#' @keywords internal
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






