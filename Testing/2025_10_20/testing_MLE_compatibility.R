source("../R/rhSSALT.R")
source("../R/MLEhSSALT.R")
source("../R/MLE_Exp.R")
source("../R/EM_algorithm_censored.R")
source("../R/EM_algorithm_interval.R")
source("../R/MLE_Geo.R")
source("../R/bs_bca_continuous_help_functions.R")
source("../R/CIbca_hSSALT.R")
source("../R/CIbs_hSSALT.R")
source("../R/CIsay_hSSALT.R")
source("../R/CIhSSALT.R")
source("../R/print.CIhSSALT.R")
source("../R/print.hSSALTMLE.R")

library(Rcpp)

sourceCpp("../src/EM_arma_int_cont.cpp")


############Variables from PDF
n <- 35
tau <- 8
theta1 <- exp(3.5)
theta21 <- exp(-0.2)
theta22 <- exp(2)
p <- 0.4
#B <- 500
r <- 30
#delta <- 0.5

data <- c(1.46,2.22,3.92,4.24,5.47,5.6,6.12,6.56,8.01,8.1,8.22,8.3,8.59,8.77,8.8,8.8,8.84,8.9,8.97,8.98,9.43,9.62,
          9.87,11.14,11.85,12.14,13.49,14.19,14.33,15.28,16.58,17.8,21.09,26.34,28.66)


resMLE.1 <- MLEhSSALT(data,n,2,r=n,tau=tau,theta21 = theta21,theta22 = theta22,p=p,language = "CPP", monitoring="continuous")
resMLE.2 <- MLEhSSALT(data,n,2,r=r,tau=tau,theta21 = theta21,theta22 = theta22,p=p,language = "CPP", monitoring="continuous")
resMLE.3 <- MLEhSSALT(data,n,1,tau=c(8,16),theta21 = theta21,theta22 = theta22,p=p,language = "CPP", monitoring="continuous")
resMLE.4 <- MLEhSSALT(data,n,1,tau=c(8,20),theta21 = theta21,theta22 = theta22,p=p,language = "CPP", monitoring="continuous")
resMLE.5 <- MLEhSSALT(data,n,1,tau=c(8,24),theta21 = theta21,theta22 = theta22,p=p,language = "CPP", monitoring="continuous")


############Given results
#Type 2 - r=n
theta1.1 <- 31.449
theta21.1 <- 0.798
theta22.1 <- 6.612
p.1 <- 0.379

abs(round(resMLE.1$mle$theta1,3) - theta1.1)
abs(round(resMLE.1$mle$theta21, 3) - theta21.1)
abs(round(resMLE.1$mle$theta22, 3) - theta22.1)
abs(round(resMLE.1$mle$p1, 3) - p.1)

#Type 2 - r=30
theta1.2 <- 31.449
theta21.2 <- 0.791
theta22.2 <- 6.447
p.2 <- 0.372

abs(round(resMLE.2$mle$theta1,3) - theta1.2)
abs(round(resMLE.2$mle$theta21, 3) - theta21.2)
abs(round(resMLE.2$mle$theta22, 3) - theta22.2)
abs(round(resMLE.2$mle$p1, 3) - p.2)

#Type 1 - tau2 = 16
theta1.3 <- 31.449
theta21.3 <- 0.818
theta22.3 <- 7.048
p.3 <- 0.395

abs(round(resMLE.3$mle$theta1,3) - theta1.3)
abs(round(resMLE.3$mle$theta21, 3) - theta21.3)
abs(round(resMLE.3$mle$theta22, 3) - theta22.3)
abs(round(resMLE.3$mle$p1, 3) - p.3)

#Type 1 - tau2=20
theta1.4 <- 31.449
theta21.4 <- 0.820
theta22.4 <- 7.085
p.4 <- 0.396

abs(round(resMLE.4$mle$theta1,3) - theta1.4)
abs(round(resMLE.4$mle$theta21, 3) - theta21.4)
abs(round(resMLE.4$mle$theta22, 3) - theta22.4)
abs(round(resMLE.4$mle$p1, 3) - p.4)

#Type 1 - tau2=24
theta1.5 <- 31.449
theta21.5 <- 0.830
theta22.5 <- 7.314
p.5 <- 0.404

abs(round(resMLE.5$mle$theta1,3) - theta1.5)
abs(round(resMLE.5$mle$theta21, 3) - theta21.5)
abs(round(resMLE.5$mle$theta22, 3) - theta22.5)
abs(round(resMLE.5$mle$p1, 3) - p.5)