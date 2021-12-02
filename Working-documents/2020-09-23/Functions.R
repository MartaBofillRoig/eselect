#####################################################################################
# Sample size calculations in terms of odds ratio  

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
# f_OR; simulates binary outcomes and computes the statistic

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
# test_f: computes the statistical test (OR)

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
# OR_function: computes the odds ratio

OR_function<- function(p0, p1){
  
  OR<- (p1/(1-p1))/(p0/(1-p0))
  return(OR)
}

##############################################################  
# fun_p0: computes the probability under the control group
# based on the pooled probability and the odds ratio

# p: pooled probability
# l: odds ratio

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


