
# check OR transformation
OR.function<- function(p0, p1){
  OR<- (p1/(1-p1))/(p0/(1-p0))
  return(OR)
}

scenarios$OR = mapply(OR.function,scenarios$p0_e1,scenarios$p1_e1) 
(scenarios$OR == scenarios$OR1)
(scenarios$OR - scenarios$OR1)