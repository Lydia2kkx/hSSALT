source("../../R/rhSSALT.R")
source("../../R/MLEhSSALT.R")
source("../../R/MLE_Exp.R")
source("../../R/EM_algorithm_censored.R")
source("../../R/EM_algorithm_interval.R")
source("../../R/MLE_Geo.R")
source("../../R/bs_bca_continuous_help_functions.R")
source("../../R/CIbca_hSSALT.R")
source("../../R/CIbs_hSSALT.R")
source("../../R/CIhSSALT.R")
source("../../R/CIsay_hSSALT.R")
source("../../R/CIhSSALT.R")
source("../../R/print.CIhSSALT.R")
source("../../R/print.hSSALTMLE.R")

library(Rcpp)
sourceCpp("../../src/EM_arma_int_cont.cpp")

# # Fixed Variables
# n <- 25
# tau <- c(8,16)
# B <- 500
# theta1 <- 30
# theta21  = 2
# theta22  = 15
# p = 0.4
# 
# language = "CPP"
# monitoring = "continuous"
# 
# set.seed(2)
# sample <- rhSSALT(n, 1, tau = tau, theta1 = theta1,theta21 = theta21, theta22 = theta22, p = p,monitoring = monitoring)
# 
# resMLE <- MLEhSSALT(sample$Censored_dat, n, 1,tau = tau, theta21 = theta21,
#                     theta22 = theta22, p = p,language = language, monitoring = monitoring)
# 
# CI <- CIhSSALT(sample$Censored_dat, n, 1, tau = tau, MLEhSSALT_Obj = resMLE, language = language,monitoring = monitoring, B = B, CImethod="asymptotic")
# CI
# 
# CI <- CIhSSALT(sample$Censored_dat, n, 1, tau = tau, MLEhSSALT_Obj = resMLE, language = language,monitoring = monitoring, B = B, CImethod="percentile")
# CI
# 
# CI <- CIhSSALT(sample$Censored_dat, n, 1, tau = tau, MLEhSSALT_Obj = resMLE, language = language,monitoring = monitoring, B = B, CImethod="bca")
# CI
# 
# 
# 
# 
# #### Type 2
# # Fixed Variables
# n <- 25
# tau <- 8
# B <- 500
# theta1 <- 30
# theta21  = 2
# theta22  = 15
# p = 0.4
# r <- 20
# 
# language = "R"
# monitoring = "continuous"
# 
# set.seed(2)
# sample <- rhSSALT(n, 2, r=r, tau = tau, theta1 = theta1,theta21 = theta21, theta22 = theta22, p = p,monitoring = monitoring)
# 
# resMLE <- MLEhSSALT(sample$Censored_dat, n, 2, r=r,tau = tau, theta21 = theta21,
#                     theta22 = theta22, p = p,language = language, monitoring = monitoring)
# 
# CI <- CIhSSALT(sample$Censored_dat, n, 2, r=r, tau = tau, MLEhSSALT_Obj = resMLE, language = language,monitoring = monitoring, B = B, CImethod="asymptotic")
# CI
# 
# CI <- CIhSSALT(sample$Censored_dat, n, 2, r=r, tau = tau, MLEhSSALT_Obj = resMLE, language = language,monitoring = monitoring, B = B, CImethod="percentile")
# CI
# 
# CI <- CIhSSALT(sample$Censored_dat, n, 2, r=r, tau = tau, MLEhSSALT_Obj = resMLE, language = language,monitoring = monitoring, B = B, CImethod="bca")
# CI




#### Interval
# Fixed Variables
n <- 35
tau <- c(8,16)
B <- 500
theta1 <- 30
theta21  = 2
theta22  = 15
p = 0.4
delta = 1

language = "R"
monitoring = "interval"

set.seed(2)
sample <- rhSSALT(n, 1, tau = tau, theta1 = theta1,theta21 = theta21, theta22 = theta22, p = p,monitoring = monitoring, delta=delta)
new_data <- sample$Full_dat

resMLE <- MLEhSSALT(new_data, n, 1,tau = tau, theta21 = theta21,
                    theta22 = theta22, p = p,language = language, monitoring = monitoring, delta=delta)

CI <- CIhSSALT(new_data, n, 1, tau = tau, MLEhSSALT_Obj = resMLE, language = language,monitoring = monitoring, B = B, CImethod="asymptotic", delta=delta)
CI

CI <- CIhSSALT(new_data, n, 1, tau = tau, MLEhSSALT_Obj = resMLE, language = language,monitoring = monitoring, B = B, CImethod="percentile", delta=delta)
CI

CI <- CIhSSALT(new_data, n, 1, tau = tau, MLEhSSALT_Obj = resMLE, language = language,monitoring = monitoring, B = B, CImethod="bca", delta=delta)
CI