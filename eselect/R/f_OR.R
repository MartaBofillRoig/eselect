##############################################################
#' f_OR
#' @description simulates binary outcomes and computes the statistic
#' @keywords internal
#' @export
#' @import stats
#' @return Two-sample test statistics using unpooled variance estimator

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
