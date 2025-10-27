source("../../R/rhSSALT.R")
source("../../R/MLEhSSALT.R")
source("../../R/MLE_Exp.R")
source("../../R/EM_algorithm_censored.R")
source("../../R/EM_algorithm_interval.R")
source("../../R/MLE_Geo.R")
source("../../R/bs_bca_continuous_help_functions.R")
source("../../R/CIbca_hSSALT.R")
source("../../R/CIbs_hSSALT.R")
source("../../R/CIsay_hSSALT.R")
source("../../R/CIhSSALT.R")
source("../../R/print.CIhSSALT.R")
source("../../R/print.hSSALTMLE.R")

source("../../R/HomohSSALT.R")
source("../../R/print.hSSALTtest.R")

library(Rcpp)

sourceCpp("../../src/EM_arma_int_cont.cpp")

n <- 35
tau <- c(8,16)
theta1 <- 30
theta21  = 2
theta22  = 15
p = 0.4
B <- 500


#Type 1
set.seed(2)
sample <- rhSSALT(n, 1, tau = tau, theta1 = theta1, theta21 = theta21, theta22 = theta22, p = p, monitoring = "continuous")
resMLE <- MLEhSSALT(sample$Censored_dat, n, 1, tau = tau, theta21 = theta21, theta22 = theta22, p = p, language = "CPP", monitoring = "continuous")

CIhSSALT(sample$Full_dat, n, 1, tau = tau, MLEhSSALT_Obj = resMLE, language = "CPP", monitoring = "continuous", B = B, CImethod = "bca")

#Type 2
set.seed(2)
sample <- rhSSALT(n, 2, r=30, tau = tau[1], theta1 = theta1, theta21 = theta21, theta22 = theta22, p = p, monitoring = "continuous")
resMLE <- MLEhSSALT(sample$Censored_dat, n, r=30, 2, tau = tau[1], theta21 = theta21, theta22 = theta22, p = p, language = "CPP", monitoring = "continuous")

CIhSSALT(sample$Censored_dat, n, r=30, 2, tau = tau[1], MLEhSSALT_Obj = resMLE, language = "CPP", monitoring = "continuous", B = B, CImethod = "bca")