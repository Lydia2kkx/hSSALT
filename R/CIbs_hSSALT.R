
CIbs_hSSALT <-function(data, n, censoring, tau , r, monitoring, delta, alpha, B, theta_21,
            theta_22, p, maxit, tol, language, parallel, ncores){
  
  same_threshold <- 0.05
  p_threshold <- 0.01
  p_adjust <- TRUE ##adjust p with same results to 0.5
  data_adjust <- TRUE ##discard same results of theta21 and theta22
  
  if (monitoring=="Interval"){
    
    
    B <- 1000
    j <- 1  ### iteration counter
    iter <- 1
    Results <- list()
    while (j <= B){
      ######First stress level
      T1_simu <- rexp(n, rate = 1/mle1)
      n1 <- sum(T1_simu <= tau[1])
      if((n1 == 0)){
        iter <- iter + 1
        next
      }else{
        T1 <- stl_sort(T1_simu[T1_simu <= tau[1]])
        n1j <- table_factor_cpp(cut(T1, breaks = seq(0, tau[1], delta)))
        theta1_hat <- uniroot(func_theta1, c(1, 1000), extendInt = "yes", n1j = n1j, n1 = n1)$root
        
        #####Second stress level
        #####Simulate the mixture
        n21 <- rbinom(1, (n-n1), p1_em)
        T21 <- stl_sort(-theta21_em*log(1-runif(n21)) + tau[1])
        n22 <- n-n1-n21
        T22 <- stl_sort(-theta22_em*log(1-runif(n22)) + tau[1])
        T2 <- stl_sort(c(T21, T22))
        n2 <- sum(T2 < tau[2])
        n2j <- table_factor_cpp(cut(T2[T2 < tau[2]], breaks = seq(tau[1], tau[2], delta)))
        data <- c(rep(j-1, n2j), rep(q2, n-n1-n2))
        d <- c(rep(1, n2), rep(0, n-n1-n2))
        q2 <- (tau[2]-tau[1])/delta
        
        parameter_starts <- data.frame(omega1=p, theta21=theta21, theta22)
        
        if(language=="CPP"){
          if(nrow(parameter_starts)==1){
            EM_mle <- EM_algorithm_interval_arma(ind =1, data = data, N = maxit, 
                                                      delta = delta, d=d, parameter_starts=parameter_starts, q2=q2, tol = tol) 
          }else if (parallel==T){
            cl <- parallel::makeCluster(getOption("cl.cores", ncores))
            EM_mle <- parLapply(cl, 1:nrow(parameter_starts), EM_algorithm_interval_arma, data = data_starts, N = maxit,
                                delta = delta, d=d, parameter_starts=parameter_starts, q2=q2, tol=tol)
            parallel::stopCluster(cl)
          }else{
            EM_mle <- lapply(1:nrow(parameter_starts), EM_algorithm_interval_arma, data = data_starts, N = maxit,
                             delta = delta, d=d, parameter_starts=parameter_starts, q2=q2, tol = tol)
          }
          
        }else{
          
          if (nrow(parameter_starts)==1){
            EM_mle <- EM_algorithm_interval(ind =1, data = data_starts, N = maxit, 
                                            delta = delta, d=d, parameter_starts=parameter_starts, q2=q2, tol = tol)
            
          }else if(parallel == TRUE){
              cl <- parallel::makeCluster(getOption("cl.cores", ncores))
              EM_mle <- parLapply(cl, 1:nrow(parameter_starts), EM_algorithm_interval, data = data_starts, N = maxit,
                                  delta = delta, d=d, parameter_starts=parameter_starts, q2=q2, tol = tol)
              parallel::stopCluster(cl)
              
            }else {
              
              EM_mle <- lapply(1:nrow(parameter_starts), EM_algorithm_interval, data = data_starts, N = maxit,
                               delta = delta, d=d, parameter_starts=parameter_starts, q2=q2, tol = tol)
            }
        #####Apply the EM algorithm
        # model_list <- lapply(1:nrow(parameter_starts), EM_algorithm_interval, data = data, delta = delta) 
        }
        Estimate_df <- do.call("rbind.data.frame", EM_mle)
        
        ###Add some variables to the result data frame
        Estimate_df$theta1 <- theta1_hat
        Estimate_df$n2 <- n2
        ### If the estimated theta22 is larger than mle1, invalid result, set prob1 (p) to NA
        Estimate_df$prob1 <- ifelse(Estimate_df$theta22 >= theta1_hat, NA, Estimate_df$prob1)
        ### If the estimated prob1 is too small (< 0.01) or too large (> 0.99), or the estimated theta21 and 
        ### theta22 are too close (less than 1% of theta22), homogeneous result 
        Estimate_df$info <- ifelse((Estimate_df$prob1 < p_threshold) | (Estimate_df$prob2 < p_threshold) | (Estimate_df$theta22 - Estimate_df$theta21 <= same_threshold*Estimate_df$theta22), 
                                   "homo", ifelse(is.na(Estimate_df$prob1), NA, "heter"))
        Estimate_df$mle2_classi <- uniroot(func_theta2, c(1, 1000), extendInt = "yes", n2j = n2j, n1 = n1, n2 = n2)$root
        if(p_adjust){
          ### if homogeneous result, set prob1 to 0.5 at first
          Estimate_df$prob1[Estimate_df$info == "homo"] <- Estimate_df$prob2[Estimate_df$info == "homo"] <- 0.5
        }
        if(data_adjust){
          ### if homogeneous result, set prob1 to NA since we only consider the heterogeneous result
          Estimate_df$prob1[Estimate_df$info == "homo"] <- NA 
        }
        if(all(is.na(Estimate_df$prob1))){ ###If all NAs in prob1, next iteration
          #seed <- seed + 1
          iter <- iter + 1
          next
        }else{ ###Otherwise, assign it to the kth sublist
          Results[[j]] <- Estimate_df
        }
      }
      #seed <- seed + 1
      iter <- iter + 1
      j <- j + 1 ### increase the iteration counter
    }
    
    ###### Omit the NA results in each sublist. For some initial points, the results might be invalid
    Results <- lapply(Results, na.omit)
    ###### Find the estimated parameters with the maximal log likelihood
    max_loglik <- lapply(Results, find_max)
    max_loglik_df <- do.call("rbind.data.frame", max_loglik)
    
    ####Find the maximum if the initial values are vectors
    #### Otherwise, it is very easy since find_max() is useless.
    #### Boostrap Percentile CIs easily take the corresponding quantiles for each parameter
    bootstrap_distri <- cbind.data.frame(i = 1:nrow(max_loglik_df), max_loglik_df)
    ### alpha is the significant level, take theta1 as an example
    # quantile(bootstrap_distri$theta1, c(alpha/2, 1-alpha/2))
    
  }else{
    
    ########## Continuous ##########
    
    if(censoring==2){
      tau <- c(tau, data[r])
    }
    
    j <- 1
    iter <- 1 ## iteration counter (how many times we need to obtain B times valid results)
    Results <- list()
    while (j <= B){
      ######First stress level
      T1_simu <- rexp(n, rate = 1/mle1)
      n1 <- sum(T1_simu <= tau[1]) # Here I show the Type-I case, just use tau[1] and tau[2] instead
      if((n1 == 0)){
        iter <- iter + 1
        next
      }else{
        T1 <- stl_sort(T1_simu[T1_simu <= tau[1]])
        theta1_hat <- sum(T1, (n-n1)*tau[1]) / n1
        n21 <- rbinom(1, (n-n1), p1_em)
        T21 <- stl_sort(-theta21_em*log(1-runif(n21)) + tau[1])
        n22 <- n-n1-n21
        T22 <- stl_sort(-theta22_em*log(1-runif(n22)) + tau[1])
        T2 <- stl_sort(c(T21, T22))
        n2 <- sum(T2 < tau[2])
        d <- c(rep(1, n2), rep(0, n-n1-n2))
        t22 <- apply(cbind(T2, tau[2]), 1, min) # observed or censored data
        
        
        
        # model_list <- lapply(1:nrow(parameter_starts), EM_algorithm_censored, data = t22 - tau[1])
        
        
        parameter_starts <- data.frame(omega1=p, theta21=theta21, theta22)
        
        
        if(language == "CPP") {
          
          if (dim(parameter_starts)[1]==1){
            EM_mle <- EM_algorithm_censored_arma(ind =1, data=t22-tau[1], d=d, N=maxit, parameter_starts = parameter_starts, tol = tol)
          
          }else if(parallel == TRUE){
            cl <- parallel::makeCluster(getOption("cl.cores", ncores))
            EM_mle <- parLapply(cl, 1:nrow(parameter_starts), EM_algorithm_censored_arma, data = t22 - tau[1], d=d, N=maxit,
                               parameter_starts = parameter_starts, tol = tol)
            parallel::stopCluster(cl)
              
          }else{
            EM_mle <- lapply(1:nrow(parameter_starts), EM_algorithm_censored_arma, data = t22 - tau[1], d=d, N=maxit,
                               parameter_starts = parameter_starts, tol = tol)
            }
            
          
        }else{
          
          if (dim(parameter_starts)[1]==1){
            EM_mle <- EM_algorithm_censored(ind =1, data=t22-tau[1], d=d, N=maxit, parameter_starts = parameter_starts, tol)
     
          }else if(parallel == TRUE){
            cl <- parallel::makeCluster(getOption("cl.cores", ncores))
            EM_mle <- parLapply(cl, 1:nrow(parameter_starts), EM_algorithm_censored, data = t22 - tau[1], d=d, N=maxit,
                               parameter_starts = parameter_starts, tol)
            parallel::stopCluster(cl)
              
          }else{
            EM_mle <- lapply(1:nrow(parameter_starts), EM_algorithm_censored, data = t22 - tau[1], d=d, N=maxit,
                               parameter_starts = parameter_starts, tol)
          }
      
        }  
          
        
        Estimate_df <- do.call("rbind.data.frame", EM_mle)
        
        
        Estimate_df$theta1 <- theta1_hat
        Estimate_df$n2 <- n2
        Estimate_df$pi1 <- ifelse(Estimate_df$theta22 >= theta1_hat, NA, Estimate_df$pi1)
        Estimate_df$info <- ifelse((Estimate_df$pi1 < p_threshold) | (Estimate_df$pi2 < p_threshold) | (Estimate_df$theta22 - Estimate_df$theta21 <= same_threshold*Estimate_df$theta22), 
                                   "homo", ifelse(is.na(Estimate_df$pi1), NA, "heter"))
        Estimate_df$mle2_classi <- (sum(T2[1:n2] - tau[1]) + (n - n1 - n2) * (tau[2] - tau[1]))/ n2
        if(p_adjust){
          Estimate_df$pi1[Estimate_df$info == "homo"] <- Estimate_df$pi2[Estimate_df$info == "homo"] <- 0.5
        }
        if(data_adjust){
          Estimate_df$pi1[Estimate_df$info == "homo"] <- NA 
        }
        if(all(is.na(Estimate_df$pi1))){
          #seed <- seed + 1
          iter <- iter + 1
          next
        }else{
          Results[[j]] <- Estimate_df  ##Also provide iter as additional information in the result
        }
      }
      #seed <- seed + 1
      iter <- iter + 1
      j <- j + 1 ### NN # increase the iteration counter
    }
    
    ####Omit NAs in the result list
    Results_omit <- lapply(Results, na.omit)
    
    ####Find the maximum if the initial values are vectors
    # find_max <- function(x){
    #   index <- which.max(x$loglik)
    #   y <- x[index, ]
    #   row.names(y) <- NULL
    #   return(y)
    # }
    
    max_loglik <- lapply(Results_omit, find_max)
    max_loglik_df <- do.call("rbind.data.frame", max_loglik)
    
    #### Otherwise, it is very easy since find_max() is useless.
    #### Boostrap Percentile CIs easily take the corresponding quantiles for each parameter
    bootstrap_distri <- cbind.data.frame(i = 1:nrow(max_loglik_df), max_loglik_df)
    ### alpha is the significant level, take theta1 as an example
    # quantile(bootstrap_distri$theta1, c(alpha/2, 1-alpha/2))
    
    
  }
  
  return(as.list(CI_theta1 = quantile(bootstrap_distri$theta1, c(alpha/2, 1-alpha/2)), 
               CI_theta21 = quantile(bootstrap_distri$theta21, c(alpha/2, 1-alpha/2)),
               CI_theta22 = quantile(bootstrap_distri$theta22, c(alpha/2, 1-alpha/2)),
               CI_p = quantile(bootstrap_distri$p, c(alpha/2, 1-alpha/2))
))
  
}
  