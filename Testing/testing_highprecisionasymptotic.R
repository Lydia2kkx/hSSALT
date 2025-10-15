source("R/rhSSALT.R")
source("R/MLEhSSALT.R")
source("R/MLE_Exp.R")
source("R/EM_algorithm_censored.R")
source("R/EM_algorithm_interval.R")
source("R/MLE_Geo.R")
source("R/bs_bca_continuous_help_functions.R")
source("R/CIbca_hSSALT.R")
source("R/CIbs_hSSALT.R")
source("R/CIhSSALT.R")
source("R/CIsay_hSSALT.R")
source("R/CIhSSALT.R")
source("R/print.CIhSSALT.R")
source("R/print.hSSALTMLE.R")

library(Rcpp)
sourceCpp("src/EM_arma_int_cont.cpp")


############Variables
n <- 20
tau=c(8,16)
theta1 <- 33.12
theta21 <- 2
theta22 <- 15
p <- 0.4
B <- 500
r <- 30
delta <- 0.5
seed <- 1:20

set.seed(2)
sample <- rhSSALT(n,1,tau=c(8,16),theta1=theta1,theta21 = theta21,theta22 = theta22,p=p, monitoring="continuous")
resMLE_cpp <- MLEhSSALT(sample$Censored_dat,n,1,tau=c(8,16),theta21 = theta21,theta22 = theta22,p=p,language = "CPP", monitoring="continuous")
CIhSSALT(sample$Censored_dat,n,1,tau=c(8,16),MLEhSSALT_Obj=resMLE_cpp,language = "CPP",monitoring="continuous", CImethod="asymptotic")

data <- sample$Censored_dat
###########Normal computation
T1 <- data[data<tau[1]]
n1 <- length(T1)
T2 <- data[data>=tau[1]]
mle1 <- (sum(T1) + tau[1]*(n-n1))/n1
mle1
p1 <- 1 - exp(-tau[1]/mle1)
p2 <- exp(-tau[1]/mle1) * (1 - p1*exp(-(tau[2]-tau[1])/theta21)  - (1-p1)*exp(-(tau[2]-tau[1])/theta22))
p3 <- 1 - p1 - p2
Cn <- 1/(1 - (1-p1)^n - (1-p2)^n + p3^n)

# add p1, p3
bias_part <- function(n){
  y <- 0
  for (i in 1:(n-1)) {
    for (k in 0:i) {
      C_ik <- (-1)^k * choose(n, i) * choose(i, k) * ((1 - p1)^(n-i) - p3^(n-i)) * (1-p1)^k
      tau_ik <- tau[1] / i * (n - i + k)
      re_ik <- C_ik *tau_ik
      y <- y + re_ik
    }
  }
  return(y)
}

V1 <- mle1^2/n1
alpha <- 0.05
theta1_approxCI_low <- mle1 - Cn * bias_part(n) - qnorm(1-alpha/2) * sqrt(V1)
theta1_approxCI_low_normal <- ifelse(theta1_approxCI_low < 0, 0L, theta1_approxCI_low)
theta1_approxCI_up_normal <- mle1 - Cn * bias_part(n) + qnorm(1-alpha/2) * sqrt(V1)
####If with high precision
T1_simu <- Rmpfr::mpfr(1, 256)*data
n1 <- sum(T1_simu <= tau[1])
T1 <- sort(T1_simu[T1_simu <= tau[1]])
mle1 <- Rmpfr::mpfr(sum(T1, (n-n1)*tau[1]) / n1, 256)

p1 <- 1 - exp(-tau[1]/mle1)
p2 <- exp(-tau[1]/mle1) * (1 - p1*exp(-(tau[2]-tau[1])/theta21) - (1-p1)*exp(-(tau[2]-tau[1])/theta22))
p3 <- 1 - p1 - p2
Cn <- 1/(1 - (1-p1)^n - (1-p2)^n + p3^n)

# add mle1, tau
bias_part_Mpfr <- function(n){
  y <- 0
  
  for (i in 1:(n-1)) {
    for (k in 0:i) {
      n <- Rmpfr::mpfr(n, 256)
      C_ik <- (-1)^k * Rmpfr::chooseMpfr(n, i) * Rmpfr::chooseMpfr(i, k) * ((1 - p1)^(n-i) - p3^(n-i)) * (1-p1)^k
      tau_ik <- tau[1] / i * (n - i + k)
      re_ik <- C_ik *tau_ik
      y <- y + re_ik
    }
  }
  return(y)
}

V1 <- mle1^2/n1
# V1

####################################################################################
######Censored
bias <- bias_part_Mpfr(n)
theta1_approxCI_low <- mle1 - Cn * bias - qnorm(1-alpha/2)*sqrt(V1)
theta1_approxCI_low <- ifelse(theta1_approxCI_low < 0, 0L, theta1_approxCI_low)
#Yao: please check if this asNumeric function can be replaced by any other basic function
theta1_approxCI_low_high <- sapply(theta1_approxCI_low[1], Rmpfr::asNumeric) 
theta1_approxCI_up <- mle1 - Cn * bias + qnorm(1-alpha/2)*sqrt(V1)
theta1_approxCI_low_high <- sapply(theta1_approxCI_up[1], Rmpfr::asNumeric)

theta1_approxCI_low_high <- Rmpfr::asNumeric(theta1_approxCI_low[1])
theta1_approxCI_up_high  <- Rmpfr::asNumeric(theta1_approxCI_up[1])


theta1_approxCI_low_normal
theta1_approxCI_low_high
theta1_approxCI_up_normal
theta1_approxCI_up_high

resMLE_cpp <- MLEhSSALT(sample$Censored_dat,n,1,tau=c(8,16),theta21 = theta21,theta22 = theta22,p=p,language = "CPP", monitoring="continuous")
CI_cpp <- CIhSSALT(sample$Censored_dat,n,1,tau=c(8,16),MLEhSSALT_Obj=resMLE_cpp,language = "CPP",monitoring="continuous", B=B, CImethod="asymptotic")
CI_cpp

set.seed(2)
sample <- rhSSALT(n,2,tau=10, r = r,theta1=theta1,theta21 = theta21,theta22 = theta22,p=p, monitoring="continuous")
MLEhSSALT(sample$Censored_dat,n,2,tau=c(8,16),theta21 = theta21,theta22 = theta22,p=p,language = "CPP", monitoring="continuous")
MLEhSSALT(sample$Censored_dat,n,2,tau=8,theta21 = theta21,theta22 = theta22,p=p,language = "CPP", monitoring="continuous")
MLEhSSALT(sample$Censored_dat,n,censoring =2,tau=10, r=r, theta21 = theta21,theta22 = theta22,p=p,language = "CPP", monitoring="continuous")
MLEhSSALT(sample$Censored_dat,n,censoring =2,tau=10, r=r, theta21 = theta21,theta22 = theta22,p=p,language = "R", monitoring="continuous")
MLEhSSALT(sample$Censored_dat,n,censoring =1,tau=c(8,16), r=r, theta21 = theta21,theta22 = theta22,p=p,language = "R", monitoring="continuous")
