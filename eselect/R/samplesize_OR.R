##############################################################
#' samplesize_OR
#' @description Sample size calculations in terms of odds ratio
#' @keywords internal
#' @export
#' @return Sample size for odds ratios

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
