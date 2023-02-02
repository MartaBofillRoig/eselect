##############################################################
#' fun_p0
#' @description computes the probability under the control group
#' based on the pooled probability and the odds ratio
#' p: pooled probability
#' l: odds ratio
#' @keywords internal
#' @export
#' @return Probability under the control group based on OR and pooled probability


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
