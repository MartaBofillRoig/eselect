##############################################################
#' OR_function
#' @description computes the odds ratio
#' @keywords internal
#' @export
#' @return Odds ratio calculation based on marginal probabilities

OR_function<- function(p0, p1){

  OR<- (p1/(1-p1))/(p0/(1-p0))
  return(OR)
}
