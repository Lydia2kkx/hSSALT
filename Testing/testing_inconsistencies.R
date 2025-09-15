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
theta21 <- 2
theta22 <- 15
p <- 0.4
B <- 500
r <- 30
delta <- 0.5
seed <- 1:20

set.seed(16)
sample <- rhSSALT(n,1,tau=tau,theta1=theta1,theta21 = theta21,theta22 = theta22,p=p, monitoring="continuous", delta=delta)

resMLE_r <- MLEhSSALT(sample$Censored_dat,n,1,tau=tau,theta21 = theta21,theta22 = theta22,p=p,language = "R", monitoring="continuous",delta=delta)
CI_r <- CIhSSALT(sample$Censored_dat,n,1,tau=tau,MLEhSSALT_Obj=resMLE_r,language = "R",monitoring = "continuous", B=B,grid=F, CImethod="asymptotic",delta=delta)

set.seed(16)
sample <- rhSSALT(n,1,tau=tau,theta1=theta1,theta21 = theta21,theta22 = theta22,p=p, monitoring="continuous",delta=delta)
resMLE_cpp <- MLEhSSALT(sample$Censored_dat,n,1,tau=tau,theta21 = theta21,theta22 = theta22,p=p,language = "CPP", monitoring="continuous",delta=delta)
CI_cpp <- CIhSSALT(sample$Censored_dat,n,1,tau=tau,MLEhSSALT_Obj=resMLE_cpp,language = "CPP",monitoring = "continuous", B=B,grid=F, CImethod="asymptotic",delta=delta)

# resMLE_r$mle$theta1
# resMLE_cpp$mle$theta1
# 
# sprintf("%.30f", resMLE_r$mle$p1)
# sprintf("%.30f", resMLE_cpp$mle$p1)

# mydfR <- read.csv("theta22R.csv")
# mydfCPP <- read.csv("theta22CPP.csv")
# mydfR$theta22
# mydfCPP$theta22

# mydfR$theta22 - mydfCPP$theta22[1:807]

CI_r
CI_cpp