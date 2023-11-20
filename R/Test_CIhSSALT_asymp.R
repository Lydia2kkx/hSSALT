library(EMcensoring)

#################
#
### Test of CIhSSalt
#
##################

# equal results

#
# censoring type I
#

tau = c(5,15)
test7 <- rhSSALT(n=n, censoring = 1, tau= tau, monitoring = "continuous", theta1 = exp(3.5), theta21 = exp(4), theta22 = exp(5), p=0.5)

output_cont_full_i = MLEhSSALT(data = test7$Full_dat, n=n, censoring = 1, tau = tau, monitoring = "continuous", theta21 = exp(1.9), 
                               theta22 = exp(3.7), p=0.7)

output_cont_cen_i = MLEhSSALT(data = test7$Censored_dat, n=n, censoring = 1, tau = tau, monitoring = "continuous", theta21 = exp(1.9), 
                              theta22 = exp(3.7), p=0.7)
cbind(output_cont_full_i,
      output_cont_cen_i)




# type I

test_CI_c_1 <- CIhSSALT(data = test7$Full_dat, censoring = 1 , n=n, tau=tau, monitoring = "continuous" ,
                    delta = NULL, CImethod = "asymptotic" , alpha = 0.05 , B = 1000 , theta1=output_cont_full_i$mle1,
                    theta21=output_cont_full_i$mle2.theta21 , theta22=output_cont_full_i$mle2.theta22 , 
                    p=output_cont_full_i$mle2.p1 , maxit=1000 , language = "CPP" , parallel = FALSE)


# singular matrix, bad


# type II

### MLE_exp Test

censoring=2
tau = c(5,15)

n=50

# test1 <- rhSSALT(n=n, censoring = censoring, tau= tau, monitoring = "continuous",delta = delta, theta1 = exp(3), theta21 = exp(1.9), theta22 = exp(3.7), p=0.7)

tau = 5
r=40
test1 <- rhSSALT(n=n, censoring = censoring, tau= tau, r=r, monitoring = "continuous", theta1 = exp(3), theta21 = exp(1.9), theta22 = exp(3.7), p=0.7)

test1$Censored_dat

test1$Full_dat

length(test1$Full_dat)


output_cont_full_ii = MLEhSSALT(data = test1$Full_dat, n=n, censoring = 2, tau = tau, r=r , monitoring = "continuous", theta21 = 4, 
                                theta22 = 7, p=0.9)

output_cont_cen_ii = MLEhSSALT(data = test1$Censored_dat, n=n, censoring = 2, tau = tau, r=r, monitoring = "continuous", theta21 = 4, 
                               theta22 = 7, p=0.9)
cbind(output_cont_full_ii,
      output_cont_cen_ii)


test_CI_c_2 <- CIhSSALT(data = test7$Full_dat, censoring = censoring , n=n, tau=tau, r=20,  monitoring = "continuous" ,
                        delta = NULL, CImethod = "asymptotic" , alpha = 0.05 , B = 1000 , theta1=output_cont_full_ii$mle1,
                        theta21=output_cont_full_ii$mle2.theta21 , theta22=output_cont_full_ii$mle2.theta22 , 
                        p=output_cont_full_i$mle2.p1 , maxit=1000 , language = "CPP" , parallel = FALSE)

# 1: In sqrt(Vp) : NaNs wurden erzeugt
# 2: In sqrt(Vp) : NaNs wurden erzeugt


###############
#
## Interval 
#
################


delta = 0.5
tau= c(2,3)
n=50

test5 <- rhSSALT(n=n, censoring = 1, tau= tau, monitoring = "interval",delta = delta, theta1 = exp(3), theta21 = exp(1.4), theta22 = exp(0.5), p=0.5)

test5$Censored_dat

test5$Full_dat

sum(test5$Full_dat) #now fifty


output_full <- MLEhSSALT(data = test5$Full_dat, n=n, censoring = 1, tau = tau, monitoring = "interval", delta = delta, theta21 = 1.5, 
                         theta22 = 3.7, p=0.4)

output_censored <- MLEhSSALT(data = test5$Censored_dat, n=n, censoring = 1, tau = tau, monitoring = "interval", delta = delta, theta21 = 1.5, 
                             theta22 = 3.7, p=0.4)


test_CI_i_1 <- CIhSSALT(data = test5$Full_dat, censoring = 1 , n=n, tau=tau, monitoring = "interval" ,
                        delta = delta, CImethod = "asymptotic" , alpha = 0.05 , B = 1000 , theta1=output_full$mle1,
                        theta21=output_full$mle2.theta21 , theta22=output_full$mle2.theta22 , 
                        p=output_full$mle2.p1 , maxit=1000 , language = "CPP" , parallel = FALSE)
# Warnmeldungen:
# 1: In sqrt(V_theta21) : NaNs wurden erzeugt
# 2: In sqrt(V_theta21) : NaNs wurden erzeugt
