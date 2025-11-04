library(foreach)
library(doParallel)

source("rhSSALT.R")
source("HomohSSALT.R")

###Fixed Variables
n <- 35
theta1 <- 30
tau <- c(8,16)
theta21  = 2
theta22  = 15
p = 0.4

TEST_LOOPS <- 100000
cl <- 16

registerDoParallel(cl)

loop_indices <- split(1:TEST_LOOPS, 1:cl)

Homogeneity_Test <- function(indices) {
  results <- list()
  
  set.seed(600)
  for (k in seq_along(indices)) {
    sample <- rhSSALT(n, 1, tau = tau, theta1 = theta1, theta21 = theta21, theta22 = theta22, p = p, monitoring = "continuous")
    try_output <- tryCatch({
      test_results <- HomohSSALT(data = sample$Censored_dat, n = n, tau = tau)
      c(
        statistic = as.numeric(test_results$statistic),
        critical_value = as.numeric(test_results$parameter["Critical Value"]),
        decision = test_results$decision,
        closeness_warning = test_results$closeness_warning
      )
    }, error = function(e) {
      c(statistic = NA, critical_value = NA, decision = NA, closeness_warning = NA)
    })
    
    results[[k]] <- try_output
  }
  results
}

result_list <- foreach(indices = loop_indices) %dopar% Homogeneity_Test(indices)

results <- do.call(c, result_list)
results_df <- as.data.frame(do.call(rbind, results))
results_df
save(results_df, file = "results_df_homo_testing.RData")