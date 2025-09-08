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

############Variables
n <- 40
tau <- c(8,16)
theta1 <- 33.12
theta21 <- c(2,3,4)
theta22 <- c(15,16,17)
p <- c(0.3,0.4,0.5)
B <- 500
r <- 30
delta <- 0.5
seed <- 1:100


set.seed(10)
sample <- rhSSALT(n,1,tau=tau,theta1=theta1,theta21 = theta21,theta22 = theta22,p=p, monitoring="continuous")

resMLE_cpp <- MLEhSSALT(sample$Censored_dat,n,1,tau=tau,theta21 = theta21,theta22 = theta22,p=p,language = "CPP", monitoring="continuous")
resMLE_cpp
resMLE_r <- MLEhSSALT(sample$Censored_dat,n,1,tau=tau,theta21 = theta21,theta22 = theta22,p=p,language = "R", monitoring="continuous")
resMLE_r