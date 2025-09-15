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
library(Rmpfr)

sourceCpp("../src/EM_arma_int_cont.cpp")

############Variables
n <- 40
tau <- c(8,16)
theta1 <- 33.12
theta21 <- 2
theta22 <- 15
p <- 0.4
B <- 500
r <- 30
delta <- 0.5
seed <- 1:100


############Testing

#for (i in seed) {
#  set.seed(i)
  set.seed(57)
#  cat("i: ", i, "\n")
  sample <- rhSSALT(n,1,tau=tau,theta1=theta1,theta21 = theta21,theta22 = theta22,p=p, monitoring="interval",delta=delta)
  resMLE <- MLEhSSALT(sample$Censored_dat,n,1,tau=tau,theta21 = theta21,theta22 = theta22,p=p,language = "CPP", monitoring="interval", delta=delta)
  resMLE.test <- MLEhSSALT(sample$Censored_dat,n,1,tau=tau,theta21 = theta21,theta22 = theta22,p=p,language = "R", monitoring="interval", delta=delta)
  resMLE
  resMLE.test
  #sprintf("%.30f", resMLE$mle$theta21)
  #sprintf("%.30f", resMLE.test$mle$theta21)
#}