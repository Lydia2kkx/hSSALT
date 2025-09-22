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
theta21 <- c(1,2,3,4)
theta22 <- c(14,15,16,17)
p <- c(0.2,0.3,0.4,0.5)

theta21 <- seq(from = 1, to = 4, length.out = 4)
theta22 <- seq(from = 14, to = 17, length.out = 4)
p <- seq(from = 0.2, to = 0.5, length.out = 4)

B <- 500
r <- 30
delta <- 0.5
seed <- 1:20


parallel::detectCores()

set.seed(10)
sample <- rhSSALT(n,1,tau=tau,theta1=theta1,theta21 = theta21,theta22 = theta22,p=p, monitoring="interval", delta=delta)
resMLE <- MLEhSSALT(sample$Censored_dat,n,1,tau=tau,theta21 = theta21,theta22 = theta22,p=p,language = "R", monitoring="interval", parallel=TRUE, ncores=4, delta=delta)

CI_r <- CIhSSALT(sample$Censored_dat,n,1,tau=tau,MLEhSSALT_Obj=resMLE,language = "R",monitoring = "interval", B=B,grid=F, CImethod="percentile",delta=delta, parallel=T)
CI_r2 <- CIhSSALT(sample$Censored_dat,n,1,tau=tau,MLEhSSALT_Obj=resMLE,language = "R",monitoring = "interval", B=B,grid=F, CImethod="percentile",delta=delta)

#system.time({resMLE_r <- MLEhSSALT(sample$Censored_dat,n,1,tau=tau,theta21 = theta21,theta22 = theta22,p=p,language = "R", monitoring="continuous", parallel=TRUE, ncores=4)})
#system.time({resMLE_r2 <- MLEhSSALT(sample$Censored_dat,n,1,tau=tau,theta21 = theta21,theta22 = theta22,p=p,language = "R", monitoring="continuous", parallel=F)})


#CI_r <- CIhSSALT(sample$Censored_dat,n,1,tau=tau,MLEhSSALT_Obj=resMLE_r,language = "R",monitoring = "continuous", B=B,grid=F, CImethod="bca",delta=delta)