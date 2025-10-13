library(foreach)
library(doParallel)

## Fixed Variables
n <- 40
tau <- c(8,16)
B <- 500
theta1 <- 30

# Testing variables
theta21  = c(2, 1.5)
theta22  = c(15, 12)
p = c(0.4, 0.5)

parameter_starts <- expand.grid(theta21, theta22, p)
colnames(parameter_starts) <- c("theta21", "theta22", "p")

times <- 2
seeds <- seq(from = 2000, by = 20000, length.out = times)

#Function:
CI_TypeI <- function(ind){
  library(hSSALT)
  theta21 <- parameter_starts[ind, "theta21"]
  theta22 <- parameter_starts[ind, "theta22"]
  p <- parameter_starts[ind, "p"]
  
  results.cont.asymptotic <- list()
  results.cont.percentile <- list()
  results.cont.bca        <- list()
  
  MLE_NotNA.cont    <- rep(0, times)
  MLE_converge.cont <- rep(0, times)
  MLE_heter.cont    <- rep(0, times)
  
  for (i in seq_along(seeds)) {
    seed_i <- seeds[i]
    
    set.seed(seed_i)
    sample <- rhSSALT(n, 1, tau = tau, theta1 = theta1, theta21 = theta21, theta22 = theta22, p = p, monitoring = "continuous")
    resMLE <- MLEhSSALT(sample$Censored_dat, n, 1, tau = tau, theta21 = theta21, theta22 = theta22, p = p,
                        language = "CPP", monitoring = "continuous")
    
    if (!is.na(resMLE$info)) {
      MLE_NotNA.cont[i]    <- TRUE
      MLE_heter.cont[i]    <- (resMLE$info == "heterogeneous")
      MLE_converge.cont[i] <- (resMLE$message == "convergent")
    } else {
      MLE_NotNA.cont[i]    <- FALSE
      MLE_converge.cont[i] <- FALSE
      MLE_heter.cont[i]    <- FALSE
    }
    
    if (is.na(resMLE$info) || resMLE$info != "heterogeneous" || resMLE$message != "convergent") {
      results.cont.asymptotic[[i]] <- list(theta1=c(0,0), theta21=c(0,0), theta22=c(0,0), p=c(0,0), B=0, j=0, alpha=0, type=0)
      results.cont.percentile[[i]] <- list(theta1=c(0,0), theta21=c(0,0), theta22=c(0,0), p=c(0,0), B=0, j=0, alpha=0, type=0)
      results.cont.bca[[i]]        <- list(theta1=c(0,0), theta21=c(0,0), theta22=c(0,0), p=c(0,0), B=0, j=0, alpha=0, type=0)
    } else {
      ## asymptotic
      set.seed(seed_i)
      sample <- rhSSALT(n, 1, tau = tau, theta1 = theta1, theta21 = theta21, theta22 = theta22, p = p, monitoring = "continuous")
      resMLE <- MLEhSSALT(sample$Censored_dat, n, 1, tau = tau, theta21 = theta21, theta22 = theta22, p = p,
                          language = "CPP", monitoring = "continuous")
      results.cont.asymptotic[[i]] <- tryCatch(
        CIhSSALT(sample$Censored_dat, n, 1, tau = tau, MLEhSSALT_Obj = resMLE, language = "CPP",
                 monitoring = "continuous", B = B, grid = FALSE),
        error = function(e) list(theta1=c(NA,NA), theta21=c(NA,NA), theta22=c(NA,NA), p=c(NA,NA), B=NA, j=NA, alpha=NA, type=NA)
      )
      
      ## percentile
      set.seed(seed_i)
      sample <- rhSSALT(n, 1, tau = tau, theta1 = theta1, theta21 = theta21, theta22 = theta22, p = p, monitoring = "continuous")
      resMLE <- MLEhSSALT(sample$Censored_dat, n, 1, tau = tau, theta21 = theta21, theta22 = theta22, p = p,
                          language = "CPP", monitoring = "continuous")
      results.cont.percentile[[i]] <- tryCatch(
        CIhSSALT(sample$Censored_dat, n, 1, tau = tau, MLEhSSALT_Obj = resMLE, language = "CPP",
                 monitoring = "continuous", B = B, CImethod = "percentile", grid = FALSE),
        error = function(e) list(theta1=c(NA,NA), theta21=c(NA,NA), theta22=c(NA,NA), p=c(NA,NA), B=NA, j=NA, alpha=NA, type=NA)
      )
      
      ## bca
      set.seed(seed_i)
      sample <- rhSSALT(n, 1, tau = tau, theta1 = theta1, theta21 = theta21, theta22 = theta22, p = p, monitoring = "continuous")
      resMLE <- MLEhSSALT(sample$Censored_dat, n, 1, tau = tau, theta21 = theta21, theta22 = theta22, p = p,
                          language = "CPP", monitoring = "continuous")
      results.cont.bca[[i]] <- tryCatch(
        CIhSSALT(sample$Censored_dat, n, 1, tau = tau, MLEhSSALT_Obj = resMLE, language = "CPP",
                 monitoring = "continuous", B = B, CImethod = "bca", grid = FALSE),
        error = function(e) list(theta1=c(NA,NA), theta21=c(NA,NA), theta22=c(NA,NA), p=c(NA,NA), B=NA, j=NA, alpha=NA, type=NA)
      )
    }
  }
  
  tmp_a <- lapply(results.cont.asymptotic, function(x) data.frame(
    theta1_l=x$theta1[1], theta1_u=x$theta1[2],
    theta21_l=x$theta21[1], theta21_u=x$theta21[2],
    theta22_l=x$theta22[1], theta22_u=x$theta22[2],
    p_l=x$p[1], p_u=x$p[2], B=x$B, j=x$j, alpha=x$alpha, type=x$type
  ))
  df.cont.asymptotic.cpp <- cbind(seed=seeds, do.call(rbind, tmp_a), MLE_NotNA.cont, MLE_converge.cont, MLE_heter.cont)
  
  tmp_p <- lapply(results.cont.percentile, function(x) data.frame(
    theta1_l=x$theta1[1], theta1_u=x$theta1[2],
    theta21_l=x$theta21[1], theta21_u=x$theta21[2],
    theta22_l=x$theta22[1], theta22_u=x$theta22[2],
    p_l=x$p[1], p_u=x$p[2], B=x$B, j=x$j, alpha=x$alpha, type=x$type
  ))
  df.cont.percentile.cpp <- cbind(seed=seeds, do.call(rbind, tmp_p), MLE_NotNA.cont, MLE_converge.cont, MLE_heter.cont)
  
  tmp_b <- lapply(results.cont.bca, function(x) data.frame(
    theta1_l=x$theta1[1], theta1_u=x$theta1[2],
    theta21_l=x$theta21[1], theta21_u=x$theta21[2],
    theta22_l=x$theta22[1], theta22_u=x$theta22[2],
    p_l=x$p[1], p_u=x$p[2], B=x$B, j=x$j, alpha=x$alpha, type=x$type
  ))
  df.cont.bca.cpp <- cbind(seed=seeds, do.call(rbind, tmp_b), MLE_NotNA.cont, MLE_converge.cont, MLE_heter.cont)
  
  big_df <- rbind(
    cbind(method="asymptotic", df.cont.asymptotic.cpp),
    cbind(method="percentile", df.cont.percentile.cpp),
    cbind(method="bca",        df.cont.bca.cpp)
  )
  big_df$theta1_start  <- theta1
  big_df$theta21_start <- theta21
  big_df$theta22_start <- theta22
  big_df$p_start <- p
  
  dir.create("csv_dir", showWarnings = FALSE)
  write.csv(big_df, file = paste0("csv_dir/TypeI_cont_set", ind, ".csv"), row.names = FALSE)
  
  return(list(params = c(theta21 = theta21, theta22 = theta22, p = p),
    df_cont_asymptotic_cpp = df.cont.asymptotic.cpp,
    df_cont_percentile_cpp = df.cont.percentile.cpp,
    df_cont_bca_cpp = df.cont.bca.cpp,
    big_df = big_df
  ))
}

parallel.function <- function(i){
  result <- CI_TypeI(i)
}

## Cores
cl <- 4
registerDoParallel(cl)

result_list <- foreach(i=1:nrow(parameter_starts)) %dopar% parallel.function(i)

all_big_df <- do.call(rbind, lapply(result_list, function(z) z$big_df))
save(result_list, all_big_df, file = "TypeI_cont_CPP_results.RData")