source("../R/rhSSALT.R")
source("../R/MLEhSSALT.R")
source("../R/MLE_Exp.R")
source("../R/EM_algorithm_censored.R")
source("../R/EM_algorithm_interval.R")
source("../R/MLE_Geo.R")
source("../R/bs_bca_continuous_help_functions.R")
source("../R/CIbca_hSSALT.R")
source("../R/CIbs_hSSALT.R")
source("../R/CIhSSALT.R")
source("../R/CIsay_hSSALT.R")
source("../R/CIhSSALT.R")
source("../R/print.CIhSSALT.R")
source("../R/print.hSSALTMLE.R")

library(Rcpp)
sourceCpp("../src/EM_arma_int_cont.cpp")

# Fixed Variables
tau <- c(8,16)
B <- 500
theta1 <- 30
theta21  = 2
theta22  = 15
p = 0.4

# Testing variables
ns <- c(20,25,30,35,40,45,50)

times <- 100
seeds <- seq(from = 1, by = 1, length.out = times)

#Function:
Theta1_testing <- function(ind){
  n <- ns[ind]
  results.cont.asymptotic <- list()
  
  for (i in seq_along(seeds)) {
    seed_i <- seeds[i]
    cat("seed:", seed_i, "\n")
    
    set.seed(seed_i)
    sample <- rhSSALT(n, 1, tau = tau, theta1 = theta1,theta21 = theta21, theta22 = theta22, p = p,monitoring = "continuous")

    resMLE <- MLEhSSALT(sample$Censored_dat, n, 1,tau = tau, theta21 = theta21,
                        theta22 = theta22, p = p,language = "CPP", monitoring = "continuous")
    
    if (is.na(resMLE$info) || resMLE$info != "heterogeneous" || resMLE$message != "convergent") {
      results.cont.asymptotic[[i]] <- list(theta1_low_diff = 0, theta1_up_diff = 0)
      next
    }

    results.cont.asymptotic[[i]] <- tryCatch({
      
      data <- sample$Censored_dat
      theta21_mle <- resMLE$mle$theta21
      theta22_mle <- resMLE$mle$theta22
      
      ########### Normal computation
      T1 <- data[data < tau[1]]
      n1 <- length(T1)
      mle1 <- (sum(T1) + tau[1] * (n - n1)) / n1
      
      p1 <- 1 - exp(-tau[1] / mle1)
      p2 <- exp(-tau[1] / mle1) * (1 - p1 * exp(-(tau[2] - tau[1]) / theta21_mle) - (1 - p1) * exp(-(tau[2] - tau[1]) / theta22_mle))
      p3 <- 1 - p1 - p2
      Cn <- 1 / (1 - (1 - p1)^n - (1 - p2)^n + p3^n)
      
      bias_part <- function(n){
        y <- 0
        for (j in 1:(n - 1)) {
          for (k in 0:j) {
            C_ik <- (-1)^k * choose(n, j) * choose(j, k) * ((1 - p1)^(n - j) - p3^(n - j)) * (1 - p1)^k
            tau_ik <- tau[1] / j * (n - j + k)
            y <- y + C_ik * tau_ik
          }
        }
        y
      }
      
      V1 <- mle1^2 / n1
      alpha <- 0.05
      theta1_approxCI_low <- mle1 - Cn * bias_part(n) - qnorm(1 - alpha / 2) * sqrt(V1)
      theta1_approxCI_low_normal <- ifelse(theta1_approxCI_low < 0, 0L, theta1_approxCI_low)
      theta1_approxCI_up_normal  <- mle1 - Cn * bias_part(n) + qnorm(1 - alpha / 2) * sqrt(V1)
      
      #### High precision
      T1_simu <- Rmpfr::mpfr(1, 256) * data
      n1 <- sum(T1_simu <= tau[1])
      T1 <- sort(T1_simu[T1_simu <= tau[1]])
      mle1 <- Rmpfr::mpfr((sum(T1) + (n - n1) * tau[1]) / n1, 256)
      
      p1 <- 1 - exp(-tau[1] / mle1)
      p2 <- exp(-tau[1] / mle1) * (1 - p1 * exp(-(tau[2] - tau[1]) / theta21_mle) - (1 - p1) * exp(-(tau[2] - tau[1]) / theta22_mle))
      p3 <- 1 - p1 - p2
      Cn <- 1 / (1 - (1 - p1)^n - (1 - p2)^n + p3^n)
      
      bias_part_Mpfr <- function(n){
        y <- Rmpfr::mpfr(0, 256)
        n_mp <- Rmpfr::mpfr(n, 256)
        for (j in 1:(n - 1)) {
          for (k in 0:j) {
            C_ik <- (-1)^k * Rmpfr::chooseMpfr(n_mp, j) * Rmpfr::chooseMpfr(j, k) *
              ((1 - p1)^(n_mp - j) - p3^(n_mp - j)) * (1 - p1)^k
            tau_ik <- Rmpfr::mpfr(tau[1] / j, 256) * (n_mp - j + k)
            y <- y + C_ik * tau_ik
          }
        }
        y
      }
      
      V1 <- mle1^2 / n1
      bias <- bias_part_Mpfr(n)
      theta1_approxCI_low <- mle1 - Cn * bias - qnorm(1 - alpha / 2) * sqrt(V1)
      theta1_approxCI_low <- ifelse(theta1_approxCI_low < 0, 0L, theta1_approxCI_low)
      theta1_approxCI_low_high <- sapply(theta1_approxCI_low[1], Rmpfr::asNumeric)
      theta1_approxCI_up <- mle1 - Cn * bias + qnorm(1 - alpha / 2) * sqrt(V1)
      theta1_approxCI_up_high <- sapply(theta1_approxCI_up[1], Rmpfr::asNumeric)
      
      list(theta1_low_diff = theta1_approxCI_low_normal - theta1_approxCI_low_high, theta1_up_diff  = theta1_approxCI_up_normal  - theta1_approxCI_up_high)
    },
    error = function(e) list(theta1_low_diff = NA, theta1_up_diff = NA))
  }
  
  dir.create("high_prec_n_csvs", showWarnings = FALSE)
  df <- data.frame(seed = seeds,theta1_low_diff = sapply(results.cont.asymptotic, function(x) x[["theta1_low_diff"]]),theta1_up_diff = sapply(results.cont.asymptotic, function(x) x[["theta1_up_diff"]]))
  write.csv(df, file = paste0("high_prec_n_csvs/high_prec_", n, ".csv"), row.names = FALSE)
}

for (i in seq_along(ns)) {
  cat("n: ", ns[i], "\n")
  Theta1_testing(i)
}
