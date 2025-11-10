

bootstrap_distribution <-function(data, n, monitoring, theta1, theta21, theta22, p, censoring,
                                  tau, r, B, delta, maxit, tol, language){
  
  same_threshold <- 0.05
  p_threshold <- 0.01
  p_adjust <- TRUE ##Adjust p with same results to 0.5
  data_adjust <- TRUE ##Discard same results of theta21 and theta22
  
  j <- 1
  iter <- 1 ## Iteration counter
  Results <- list()
  
  while (j <= B){
    
    sample <- suppressWarnings(rhSSALT(n, censoring = censoring, r = r, tau = tau, theta1 = theta1, theta21 = theta21, theta22 = theta22,
                                       p = p, monitoring = monitoring, delta = delta))
    n1 <- sample$Censored_num_level[1]
    n2 <- sample$Censored_num_level[2]
    
    if((n1 == 0)){
      iter <- iter + 1
      next
    }
    
    MLE_results <- suppressWarnings(MLEhSSALT(sample$Censored_dat, n, censoring = censoring, r = r, tau = tau, theta21 = theta21, 
                            theta22 = theta22, p = p, language = language,
                            monitoring = monitoring, delta = delta))
    
    Estimate_df <- MLE_results$mle
    Estimate_df$loglik <- MLE_results$loglik
    Estimate_df$n2 <- n2
    
    if (MLE_results$message == "not convergent") {
      iter <- iter + 1
      next
    }
    
    if(is.na(Estimate_df$p1) || is.na(Estimate_df$theta1) || is.na(Estimate_df$theta22)){
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
      Results[[j]] <- Estimate_df
    }
    iter <- iter + 1
    j <- j + 1
  }
  
  ####Omit NAs in the result list
  Results_omit <- lapply(Results, na.omit)
  
  max_loglik <- lapply(Results_omit, find_max)
  max_loglik_df <- do.call("rbind.data.frame", max_loglik)
  
  return(list(cbind.data.frame(i = 1:nrow(max_loglik_df), max_loglik_df), iterations = iter - 1))
  
}
