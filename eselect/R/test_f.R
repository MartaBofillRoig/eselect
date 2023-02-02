##############################################################
#' test_f
#' @description computes the statistical test (OR)
#' @keywords internal
#' @export
#' @return Two-sample test statistics for a single endpoint

test_f <- function(OR,p0,n){

  # OR: odds ratio
  # p0: probability under control group
  # n: sample size per group

  p1 <- (OR*p0/(1-p0))/(1+(OR*p0/(1-p0)))
  den <- 1/(n*(p0*(1-p0))) + 1/(n*(p1*(1-p1)))

  test <- (log(OR))/sqrt(den)
  return(test)

}
