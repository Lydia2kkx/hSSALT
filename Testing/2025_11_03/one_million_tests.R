library(foreach)
library(doParallel)

source("rhSSALT.R")
source("HomohSSALT.R")

###Fixed Variables
n <- 35
tau <- c(8,16)
theta1 <- 30

###Grid
theta21_vals <- c(2, 4)
theta22_vals <- c(15, 20)
p_vals <- c(0.4, 0.7)
param_grid <- expand.grid(theta21 = theta21_vals,theta22 = theta22_vals,p = p_vals)

param_grid_16 <- param_grid[rep(1:nrow(param_grid), each = 2), ]

TEST_LOOPS <- 25000
dir.create("results_csv")

cl <- 16
registerDoParallel(cl)

Homogeneity_Test <- function(theta21, theta22, p) {
  results <- vector("list", TEST_LOOPS)
  for (k in seq_len(TEST_LOOPS)) {
    sample <- rhSSALT(n, 1, tau = tau, theta1 = theta1,theta21 = theta21, theta22 = theta22, p = p, monitoring = "continuous")
    try_output <- tryCatch({
      test_results <- HomohSSALT(data = sample$Censored_dat, n = n, tau = tau, M = 10000)
      c(statistic = as.numeric(test_results$statistic),
        critical_value = as.numeric(test_results$parameter["Critical Value"]),
        decision = test_results$decision,
        closeness_warning = test_results$closeness_warning)
    }, error = function(e) {
      c(statistic = NA, critical_value = NA,decision = NA, closeness_warning = NA)
    })
    results[[k]] <- try_output
  }
  results
}

result_list <- foreach(i = 1:nrow(param_grid_16)) %dopar% {
  params <- param_grid_16[i, ]
  res <- Homogeneity_Test(params$theta21, params$theta22, params$p)
  
  df <- as.data.frame(do.call(rbind, res))
  df$theta21 <- params$theta21
  df$theta22 <- params$theta22
  df$p <- params$p
  
  rep_id <- (i - 1) %% 2 + 1
  fname <- paste0("results_csv/results_theta21_", params$theta1,"_theta22_", params$theta22, "_p_", params$p, "_run_", rep_id, ".csv")
  write.csv(df, file = fname, row.names = FALSE)
  df
}

results_df <- as.data.frame(do.call(rbind, result_list))
save(results_df, file = "all_results.RData")