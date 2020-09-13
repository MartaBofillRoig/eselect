

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
  TestOR_unpooled = log((phat_group1/(1-phat_group1))/(phat_group0/(1-phat_group0)))*((1/(phat_group0*(1-phat_group0))+ 1/(phat_group1*(1-phat_group1)))/samplesize)^(-1/2)
  
  return(TestOR_unpooled)
}