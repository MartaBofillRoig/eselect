

install.packages("mvtnorm")
library(mvtnorm)

##  

alpha= 0.05
z.alpha <- qnorm(1-alpha,0,1) 

rho=0.8
corr_struct=diag(2)
corr_struct[1,2]=rho
corr_struct[2,1]=rho
pmvnorm(lower=c(-Inf,-Inf), upper=c(z.alpha,z.alpha), mean=c(0,0), corr=corr_struct)
quant <- qmvnorm(0.95, corr = corr_struct, tail = "lower.tail")
pmvnorm(lower=c(-Inf,-Inf), upper=c(quant$quantile,quant$quantile), mean=c(0,0), corr=corr_struct)


