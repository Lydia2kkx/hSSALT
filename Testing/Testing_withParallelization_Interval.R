library(foreach)
library(doParallel)

## Fixed Variables
n <- 40
tau <- c(8,16)
B <- 500
theta1 <- 30
delta <- 0.5

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
  
  results.int.asymptotic <- list()
  results.int.percentile <- list()
  results.int.bca        <- list()
  
  MLE_NotNA.int    <- rep(0, times)
  MLE_converge.int <- rep(0, times)
  MLE_heter.int    <- rep(0, times)
  
  for (i in seq_along(seeds)) {
    seed_i <- seeds[i]
    
    set.seed(seed_i)
    sample <- rhSSALT(n, 1, tau = tau, theta1 = theta1, theta21 = theta21, theta22 = theta22, p = p, monitoring = "interval", delta = delta)
    resMLE <- MLEhSSALT(sample$Censored_dat, n, 1, tau = tau, theta21 = theta21, theta22 = theta22, p = p,
                        language = "CPP", monitoring = "interval", delta = delta)
    
    if (!is.na(resMLE$info)) {
      MLE_NotNA.int[i]    <- TRUE
      MLE_heter.int[i]    <- (resMLE$info == "heterogeneous")
      MLE_converge.int[i] <- (resMLE$message == "convergent")
    } else {
      MLE_NotNA.int[i]    <- FALSE
      MLE_converge.int[i] <- FALSE
      MLE_heter.int[i]    <- FALSE
    }
    
    if (is.na(resMLE$info) || resMLE$info != "heterogeneous" || resMLE$message != "convergent") {
      results.int.asymptotic[[i]] <- list(theta1=c(0,0), theta21=c(0,0), theta22=c(0,0), p=c(0,0), B=0, j=0, alpha=0, type=0)
      results.int.percentile[[i]] <- list(theta1=c(0,0), theta21=c(0,0), theta22=c(0,0), p=c(0,0), B=0, j=0, alpha=0, type=0)
      results.int.bca[[i]]        <- list(theta1=c(0,0), theta21=c(0,0), theta22=c(0,0), p=c(0,0), B=0, j=0, alpha=0, type=0)
    } else {
      ## asymptotic
      set.seed(seed_i)
      sample <- rhSSALT(n, 1, tau = tau, theta1 = theta1, theta21 = theta21, theta22 = theta22, p = p, monitoring = "interval", delta = delta)
      resMLE <- MLEhSSALT(sample$Censored_dat, n, 1, tau = tau, theta21 = theta21, theta22 = theta22, p = p,
                          language = "CPP", monitoring = "interval", delta = delta)
      results.int.asymptotic[[i]] <- tryCatch(
        CIhSSALT(sample$Censored_dat, n, 1, tau = tau, MLEhSSALT_Obj = resMLE, language = "CPP",
                 monitoring = "interval", B = B, grid = FALSE, delta = delta),
        error = function(e) list(theta1=c(NA,NA), theta21=c(NA,NA), theta22=c(NA,NA), p=c(NA,NA), B=NA, j=NA, alpha=NA, type=NA)
      )
      
      ## percentile
      set.seed(seed_i)
      sample <- rhSSALT(n, 1, tau = tau, theta1 = theta1, theta21 = theta21, theta22 = theta22, p = p, monitoring = "interval", delta = delta)
      resMLE <- MLEhSSALT(sample$Censored_dat, n, 1, tau = tau, theta21 = theta21, theta22 = theta22, p = p,
                          language = "CPP", monitoring = "interval", delta = delta)
      results.int.percentile[[i]] <- tryCatch(
        CIhSSALT(sample$Censored_dat, n, 1, tau = tau, MLEhSSALT_Obj = resMLE, language = "CPP",
                 monitoring = "interval", B = B, CImethod = "percentile", grid = FALSE, delta = delta),
        error = function(e) list(theta1=c(NA,NA), theta21=c(NA,NA), theta22=c(NA,NA), p=c(NA,NA), B=NA, j=NA, alpha=NA, type=NA)
      )
      
      ## bca
      set.seed(seed_i)
      sample <- rhSSALT(n, 1, tau = tau, theta1 = theta1, theta21 = theta21, theta22 = theta22, p = p, monitoring = "interval", delta = delta)
      resMLE <- MLEhSSALT(sample$Censored_dat, n, 1, tau = tau, theta21 = theta21, theta22 = theta22, p = p,
                          language = "CPP", monitoring = "interval", delta = delta)
      results.int.bca[[i]] <- tryCatch(
        CIhSSALT(sample$Censored_dat, n, 1, tau = tau, MLEhSSALT_Obj = resMLE, language = "CPP",
                 monitoring = "interval", B = B, CImethod = "bca", grid = FALSE, delta = delta),
        error = function(e) list(theta1=c(NA,NA), theta21=c(NA,NA), theta22=c(NA,NA), p=c(NA,NA), B=NA, j=NA, alpha=NA, type=NA)
      )
    }
  }
  
  tmp_a <- lapply(results.int.asymptotic, function(x) data.frame(
    theta1_l=x$theta1[1], theta1_u=x$theta1[2],
    theta21_l=x$theta21[1], theta21_u=x$theta21[2],
    theta22_l=x$theta22[1], theta22_u=x$theta22[2],
    p_l=x$p[1], p_u=x$p[2], B=x$B, j=x$j, alpha=x$alpha, type=x$type
  ))
  df.int.asymptotic.cpp <- cbind(seed=seeds, do.call(rbind, tmp_a), MLE_NotNA.int, MLE_converge.int, MLE_heter.int)
  
  tmp_p <- lapply(results.int.percentile, function(x) data.frame(
    theta1_l=x$theta1[1], theta1_u=x$theta1[2],
    theta21_l=x$theta21[1], theta21_u=x$theta21[2],
    theta22_l=x$theta22[1], theta22_u=x$theta22[2],
    p_l=x$p[1], p_u=x$p[2], B=x$B, j=x$j, alpha=x$alpha, type=x$type
  ))
  df.int.percentile.cpp <- cbind(seed=seeds, do.call(rbind, tmp_p), MLE_NotNA.int, MLE_converge.int, MLE_heter.int)
  
  tmp_b <- lapply(results.int.bca, function(x) data.frame(
    theta1_l=x$theta1[1], theta1_u=x$theta1[2],
    theta21_l=x$theta21[1], theta21_u=x$theta21[2],
    theta22_l=x$theta22[1], theta22_u=x$theta22[2],
    p_l=x$p[1], p_u=x$p[2], B=x$B, j=x$j, alpha=x$alpha, type=x$type
  ))
  df.int.bca.cpp <- cbind(seed=seeds, do.call(rbind, tmp_b), MLE_NotNA.int, MLE_converge.int, MLE_heter.int)
  
  big_df <- rbind(
    cbind(method="asymptotic", df.int.asymptotic.cpp),
    cbind(method="percentile", df.int.percentile.cpp),
    cbind(method="bca",        df.int.bca.cpp)
  )
  big_df$theta1_start  <- theta1
  big_df$theta21_start <- theta21
  big_df$theta22_start <- theta22
  big_df$p_start <- p
  
  dir.create("csv_dir", showWarnings = FALSE)
  write.csv(big_df, file = paste0("csv_dir/TypeI_int_set", ind, ".csv"), row.names = FALSE)
  
  return(list(params = c(theta21 = theta21, theta22 = theta22, p = p),
    df_int_asymptotic_cpp = df.int.asymptotic.cpp,
    df_int_percentile_cpp = df.int.percentile.cpp,
    df_int_bca_cpp = df.int.bca.cpp,
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
save(result_list, all_big_df, file = "TypeI_int_CPP_results.RData")