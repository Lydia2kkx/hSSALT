

bootstrap_distribution <-function(data, n, monitoring, theta1, theta21, theta22, p, censoring,
                                  tau, r, B, delta, maxit, tol, language, parallel, ncores, grid){
  
  if(censoring == 2){
    tau <- c(tau[1], data[r])
  }
  
  #Avner: Simple Grid for now
  if (grid) {
    p_grid <- c(p*0.5,p,p*0.5+0.5)
    theta21_grid <- c(theta21*0.5,theta21,theta21*1.5)
    theta22_grid <- c(theta22*0.5,theta22,theta22*1.5)
  } else {
    p_grid <- p
    theta21_grid <- theta21
    theta22_grid <- theta22
  }
  
  same_threshold <- 0.05
  p_threshold <- 0.01
  p_adjust <- TRUE ##adjust p with same results to 0.5
  data_adjust <- TRUE ##discard same results of theta21 and theta22
  
  j <- 1
  iter <- 1 ## iteration counter (how many times we need to obtain B times valid results)
  Results <- list()
  
  while (j <= B){
    
    sample <- suppressWarnings(rhSSALT(n, 1, tau = tau, theta1 = theta1, theta21 = theta21, theta22 = theta22,
                      p = p, monitoring = monitoring, delta = delta))
    n1 <- sample$Censored_num_level[1]
    n2 <- sample$Censored_num_level[2]
    
    if((n1 == 0)){
      iter <- iter + 1
      next
    }
    
    MLE_results <- suppressWarnings(MLEhSSALT(sample$Censored_dat, n, 1, tau = tau, theta21 = theta21_grid, 
                            theta22 = theta22_grid, p = p_grid, language = language,
                            monitoring = monitoring, delta = delta, parallel = parallel, ncores = ncores))
    Estimate_df <- MLE_results$mle
    Estimate_df$loglik <- MLE_results$loglik
    Estimate_df$n2 <- n2
    
    if (MLE_results$message == "not convergent") {
      iter <- iter + 1
      next
    }
    
    if(is.na(Estimate_df$p1) || is.na(Estimate_df$theta1) || is.na(Estimate_df$theta22)){ ###If all NAs in prob1, next iteration
      iter <- iter + 1
      next
    }
    if (Estimate_df$theta22 >= Estimate_df$theta1) {
      Estimate_df$p1 <- NA
      iter <- iter + 1
      next
    }
    
    
    if ((Estimate_df$p1 < p_threshold) | (Estimate_df$p2 < p_threshold) | (Estimate_df$theta22 - Estimate_df$theta21 <= same_threshold*Estimate_df$theta22)) {
      Estimate_df$info <- "homo"
    } else if (is.na(Estimate_df$p1)) {
      Estimate_df$info <- NA
    } else {
      Estimate_df$info <- "heter"
    }
    
    if(p_adjust){
      Estimate_df$p1[Estimate_df$info == "homo"] <- Estimate_df$p2[Estimate_df$info == "homo"] <- 0.5
    }
    if(data_adjust){
      Estimate_df$p1[Estimate_df$info == "homo"] <- NA
    }
    if(all(is.na(Estimate_df$p1))){
      iter <- iter + 1
      next
    }else{
      Results[[j]] <- Estimate_df  ##Also provide iter as additional information in the result
    }
    iter <- iter + 1
    j <- j + 1 ### increase the iteration counter
    
  }
  
  ####Omit NAs in the result list
  Results_omit <- lapply(Results, na.omit)
  
  max_loglik <- lapply(Results_omit, find_max)
  max_loglik_df <- do.call("rbind.data.frame", max_loglik)
  
  #### Otherwise, it is very easy since find_max() is useless.
  #### Boostrap Percentile CIs easily take the corresponding quantiles for each parameter
  # bootstrap_distri <- cbind.data.frame(i = 1:nrow(max_loglik_df), max_loglik_df)
  
  return(list(cbind.data.frame(i = 1:nrow(max_loglik_df), max_loglik_df), iterations = iter - 1))
  
}
