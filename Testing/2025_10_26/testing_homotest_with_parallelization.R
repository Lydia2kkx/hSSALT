library(foreach)
library(doParallel)

source("../../R/rhSSALT.R")
source("../../R/HomohSSALT.R")

###Fixed Variables
n <- 35
tau <- c(8,16)
theta1 <- 30
theta21  = 2
theta22  = 15
p = 0.4
TEST_LOOPS <- 2500

set.seed(2)
sample <- rhSSALT(n, 1, tau = tau, theta1 = theta1,theta21 = theta21, theta22 = theta22,p = p, monitoring = "continuous")

cl <- 10
registerDoParallel(cl)

loop_indices <- split(1:TEST_LOOPS, 1:cl)

Homogeneity_Test <- function(indices) {
  results_10 <- numeric(length(indices))
  results_100 <- numeric(length(indices))
  
  for (k in seq_along(indices)) {
    i <- indices[k]
    
    test_10 <- HomohSSALT(data = sample$Censored_dat, n = n, tau = tau)
    val10 <- as.numeric(sprintf("%.30f", test_10$parameter[1]))
    
    test_100 <- HomohSSALT(data = sample$Censored_dat, n = n, tau = tau, M = 100000)
    val100 <- as.numeric(sprintf("%.30f", test_100$parameter[1]))
    
    results_10[k] <- val10
    results_100[k] <- val100
  }
  
  list(results_10 = results_10, results_100 = results_100)
}

result_list <- foreach(indices = loop_indices) %dopar% Homogeneity_Test(indices)

results_10 <- unlist(lapply(result_list, function(x) x[["results_10"]]))
results_100 <- unlist(lapply(result_list, function(x) x[["results_100"]]))

results_df <- data.frame(results_10, results_100)
#save(results_df, file = "results_df.RData")

results_df

# #Boxplot
# boxplot(results_df$results_10, results_df$results_100,names = c("10k", "100k"),col = c("lightgreen", "pink"))



# 
# system.time(test_10 <- HomohSSALT(data = sample$Censored_dat, n = n, tau = tau))
# system.time(HomohSSALT(data = sample$Censored_dat, n = n, tau = tau, M = 100000))