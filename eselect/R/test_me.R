##############################################################
#' test_me
#' @description computes the statistical tests (OR) for two endpoints
#' @keywords internal
#' @export
#' @return Two-sample test statistics for endpoint 1 and 2

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













