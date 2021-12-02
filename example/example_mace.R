

# Clinical outcomes at nine months:

# Probabilities in the control group and correlation

p0_e1 = 0.037 #death
p0_e2 = 0.147 #mi
p0_e3 = 0.021 #stent

# probabilities e1 and e2; and e2 and e3; prob(e1 and e3) is assumed equal 0
p0_e12 = 15/969
p0_e23 = 19/969

# 164-(36+142+20)  
  
p0_ce =  p0_e1+p0_e2+p0_e3-p0_e12-p0_e23


corr0 = ((p0_e1+p0_e2-p0_ce)-p0_e1*p0_e2)/sqrt(p0_e1*(1-p0_e1)*p0_e2*(1-p0_e2))
corr0

# 
# p0_e1 = 0.037
# p0_e2 = 0.147
# p0_e3 = 0.021
# 
# p0_ce = 0.169#0.150
# 
# 
# corr0 = ((p0_e1+p0_e2-p0_ce)-p0_e1*p0_e2)/sqrt(p0_e1*(1-p0_e1)*p0_e2*(1-p0_e2))
# corr0

# Probabilities in the treat group and correlation
# 
# p1_e1 = 0.047
# p1_e2 = 0.014
# p1_e3 = 0.035
# 
# p1_ce = 0.085
# 
# 
# corr1 = ((p1_e1+p1_e2-p1_ce)-p1_e1*p1_e2)/sqrt(p1_e1*(1-p1_e1)*p1_e2*(1-p1_e2))
# corr1


# Comparison e1 vs ce

# db <- table 2x2 results

eselect(db,p0_e1,OR1=07,p0_e2,OR2=0.8,criteria="SS",alpha=0.05,beta=0.2)

eselectsim(ss_arm=600,p0_e1=p0_e1,OR1=0.8,p0_e2=p0_e2,OR2=0.7,p0_ce=0.2,p_init=1,criteria="SS",SS_r=FALSE,alpha=0.05,beta=0.2)
