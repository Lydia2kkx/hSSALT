
CIbca_hSSALT<- function(data, n, censoring, tau, r, monitoring, delta, alpha, B,
                        theta1, theta21, theta22, p, maxit, tol, language, parallel, ncores, grid){
  
  #Avner: I'm not sure if this is the right way to go about dealing with Type-II censoring
  if(censoring==2){
    tau <- c(tau[1], data[r])
  }

  bootstrap_distri <- bootstrap_distribution(data, n, monitoring = monitoring, theta1 = theta1, theta21 = theta21,
                                             theta22 = theta22, p = p, censoring, tau, r, B, delta = delta, maxit, tol, 
                                             language, parallel, ncores, grid)

  ### alpha is the significant level, take theta1 as an example
  # quantile(bootstrap_distri$theta1, c(alpha/2, 1-alpha/2))

  ############Boostrap BCa CIs
  ### Compute the bias correction
  ### order the bootstrap estimates
  theta1_bootstrap_ordered <- sort(bootstrap_distri$theta1)
  p_bootstrap_ordered <- sort(bootstrap_distri$p1)
  theta21_bootstrap_ordered <- sort(bootstrap_distri$theta21)
  theta22_bootstrap_ordered <- sort(bootstrap_distri$theta22)

  z01 <- qnorm(sum(theta1_bootstrap_ordered < theta1) / B)
  z0p <- qnorm(sum(p_bootstrap_ordered < p) / B)
  z021 <- qnorm(sum(theta21_bootstrap_ordered < theta21) / B)
  z022 <- qnorm(sum(theta22_bootstrap_ordered < theta22) / B)

  ### BCa method is computed based on the sample data
  ### Compute the acceleration ai
  ####Calculation of a1
  ### T1 is the failures under the first stress level (simulated or real ??)
  ### Again, take Type-I as an example
  
  parameter_starts <- data.frame(omega1=p, theta21=theta21, theta22=theta22)
  
  if (monitoring == "continuous") {
    
    T1 <- data[data<tau[1]]
    n1 <- length(T1)
    T2 <- data[data>=tau[1]] 
    n2 <- length(T2)
    
    mle1ij <- c()
    for (i in 1 : n1) {
      T1_new <- T1[-i]
      mle_1new <- (sum(T1_new) + (n - n1) * tau[1]) / (n1 - 1)
      mle1ij <- c(mle1ij, mle_1new)
    }
    
    mle1ij <- c(mle1ij, rep((sum(T1) + (n - 1 - n1) * tau[1]) / n1, n2))
    mle1i <- mean(mle1ij)
    a1 <- sum((mle1i - mle1ij)^3) / (6 * (sum((mle1i - mle1ij)^2))^1.5)
    
    #########################################################################################
    # This section needs rework
    
    ####Calculation of a2s
    ### T2 is the failures under the second stress level
    Re <- list()
    for (i in 1 : n2) { #Yao: I changed (length(T2) - (n-n1-n2)) into n2, since the original one
                        #is length(T2_combine), which includes also the censored items. Thus, we
                        #need length(T2) - (n-n1-n2) for the number of observed failure under the 
                        #second stress level.
      T2_new <- T2[-i]  
      d <- c(1*(T2_new <= tau[2]), rep(0, n-n1-n2)) #Yao: I also add the censored part here
      t22_new <- apply(cbind(T2_new, tau[2]), 1, min) # observed or censored data
      if(language == "CPP"){
        if(parallel == TRUE){
          cl <- parallel::makeCluster(getOption("cl.cores", ncores))
          model_list_new <- parallel::parLapply(cl,1:nrow(parameter_starts), EM_algorithm_censored_arma, data = t22_new - tau[1], d=d, N=maxit,
                                      parameter_starts = parameter_starts, tol = tol)
          parallel::stopCluster(cl)
        }else{
          model_list_new <- lapply(1:nrow(parameter_starts), EM_algorithm_censored_arma, data = t22_new - tau[1], d=d, N=maxit,
                                   parameter_starts = parameter_starts, tol = tol)
        }
      }else{
        if(parallel == TRUE){
          cl <- parallel::makeCluster(getOption("cl.cores", ncores))
          model_list_new <- parallel::parLapply(cl,1:nrow(parameter_starts), EM_algorithm_censored, data = t22_new - tau[1], d=d, N=maxit,
                                      parameter_starts = parameter_starts, tol = tol)
          parallel::stopCluster(cl)
        }else{
          model_list_new <- lapply(1:nrow(parameter_starts), EM_algorithm_censored, data = t22_new - tau[1], d=d, N=maxit,
                                   parameter_starts = parameter_starts, tol = tol)
        }
        
      }
      
      Estimate_df_new <- data.frame(matrix(nrow = length(model_list_new), ncol = 7))
      
      posterior_l_new <- vector(mode = "list", length=length(model_list_new))
      for (j in 1:length(model_list_new)){
        Estimate_df_new[j,] <- model_list_new[[j]]$results
        posterior_l_new[[j]] <- model_list_new[[j]]$posterior
      }
      
      Re[[i]] <- Estimate_df_new
    }
  } else {
    
    ### q1 is the number of intervals under the first stress level
    q1 <- tau[1]/delta
    n1j <- data[1:q1]
    
    n1 <- sum(n1j)
    
    ### tau1j are the right bounds of intervals
    tau1j <- seq(delta, tau[1], delta)
    ### tau1j0 are the left bounds of intervals
    tau1j0 <- tau1j - delta
    
    ### tau1j are the right bounds of intervals
    tau2j <- seq(tau[1]+delta, tau[2], delta)
    ### tau1j0 are the left bounds of intervals
    tau2j0 <- tau2j - delta
    
    ### q2 is the number of intervals under the second stress level
    q2 <- (tau[2]-tau[1])/delta
    
    if (n > sum(data)) { #data is censored with a tau2 that is smaller than the input tau2
      n2j <- data[-(1:q1)] # censored
      if(q2 > length(n2j)){n2j <- c(n2j, rep(0, q2-length(n2j)))}
    }else{n2j <- data[-c((1:q1),((q1+q2+1):length(data)))]}
    
    n2 <- sum(n2j)
    
    n21 <- round(p*(n-n1))
    n22 <- sum(data)-n1-n21
    
    jackknife_samples_interval <- list()
    index <- 1
    for (bin in seq_along(data)) {
      if (data[bin] > 0) {
        for (k in 1:data[bin]) {
          temp <- data
          temp[bin] <- temp[bin] - 1
          jackknife_samples_interval[[index]] <- temp
          index <- index + 1
        }
      }
    }

    jackknife1 <- lapply(jackknife_samples_interval[1:n1], function(tbl) {
      upper_bounds <- as.numeric(sub("\\((.+),(.+)\\]", "\\2", names(tbl)))
      tbl[upper_bounds <= tau[1]]
    })
    jackknife2 <- lapply(jackknife_samples_interval[(n1+1):(n1+n2)], function(tbl) {
      upper_bounds <- as.numeric(sub("\\((.+),(.+)\\]", "\\2", names(tbl)))
      lower_bounds <- as.numeric(sub("\\((.+),(.+)\\]", "\\1", names(tbl)))
      tbl[upper_bounds <= tau[2]]
      tbl[lower_bounds >= tau[1]]
    })
    
    mle1ij <- c()
    for (i in 1 : n1) {
      n1j <- jackknife1[[i]]
      tau1j <- seq(delta, tau[1], delta)
      tau1j0 <- tau1j - delta
      mle_1new <- uniroot(func_theta1_int, c(1, 1000), extendInt = "yes", n=n,n1j = n1j, n1 = n1-1, tau1 = tau[1],tau1j0 = tau1j0, tau1j = tau1j)$root
      mle1ij <- c(mle1ij, mle_1new)
    }
    
    n1j <- data[1:q1]
    
    mle1_o <- uniroot(func_theta1_int, c(1, 1000), extendInt = "yes",n=n, n1j = n1j, n1 = n1+1, tau1j0 = tau1j0, tau1j = tau1j, tau1 = tau[1])$root
    mle1ij <- c(mle1ij, rep(mle1_o, n2))
    mle1i <- mean(mle1ij)
    a1 <- sum((mle1i - mle1ij)^3) / (6 * (sum((mle1i - mle1ij)^2))^1.5)
    
    ####Calculation of a2s
    Re <- list()
    for (i in 1 : n2) { #Yao: This (length(T2_combine) - (n-n1-n2)) should also be changed
      T2_new <- jackknife2[[i]]
      d <- c(1*(T2_new < tau[2]), rep(0, n-n1-n2))
      data_new <- c(rep((1:q2)-1, n2j), rep(q2, max(0,n-1-n1-sum(d))))
      model_list_new <- lapply(1:nrow(parameter_starts), EM_algorithm_interval, data = data_new, delta = delta, q2 = q2, d=d, N=maxit,
                               parameter_starts = parameter_starts, tol = tol)

      
      Estimate_df_new <- data.frame(matrix(nrow = length(model_list_new), ncol = 7))
      
      
      posterior_l_new <- vector(mode = "list", length=length(model_list_new))
      for (j in 1:length(model_list_new)){
        Estimate_df_new[j,] <- model_list_new[[j]]$results
        posterior_l_new[[j]] <- model_list_new[[j]]$posterior
      }
      
      Re[[i]] <- Estimate_df_new
    }
  }

  #########################################################################################

  ### Again here, if the initial values are vectors, use find_max
  ### Otherwise, simply skip
  #stab_re <- lapply(Re, check_stability)
  #stab_re_df <- do.call("rbind.data.frame", stab_re)
  #max_loglik_re <- lapply(Re, find_max)
  max_loglik_re_df <- do.call("rbind.data.frame", Re)
  
  colnames(max_loglik_re_df) <- colnames(model_list_new[[1]]$results)


  pij <- c(rep(p, n1), max_loglik_re_df$p1)
  pi <- mean(pij)
  ap <- sum((pi - pij)^3) / (6 * (sum((pi - pij)^2))^1.5)

  a21ij <- c(rep(theta21, n1), max_loglik_re_df$theta21)
  a21i <- mean(a21ij)
  a21 <- sum((a21i - a21ij)^3) / (6 * (sum((a21i - a21ij)^2))^1.5)

  a22ij <- c(rep(theta22, n1), max_loglik_re_df$theta22)
  a22i <- mean(a22ij)
  a22 <- sum((a22i - a22ij)^3) / (6 * (sum((a22i - a22ij)^2))^1.5)

  #### Compute adjusted quantiles
  ##################for theta1
  alpha11 <- pnorm(z01 + (z01 + qnorm(alpha/2)) / (1 - a1 * (z01 + qnorm(alpha/2))))
  alpha21 <- pnorm(z01 + (z01 + qnorm(1 - alpha/2)) / (1 - a1 * (z01 + qnorm(1 - alpha/2))))

  ##################for p
  alpha1p <- pnorm(z0p + (z0p + qnorm(alpha/2)) / (1 - ap * (z0p + qnorm(alpha/2))))
  alpha2p <- pnorm(z0p + (z0p + qnorm(1 - alpha/2)) / (1 - ap * (z0p + qnorm(1 - alpha/2))))

  ##################for theta2
  alpha121 <- pnorm(z021 + (z021 + qnorm(alpha/2)) / (1 - a21 * (z021 + qnorm(alpha/2))))
  alpha221 <- pnorm(z021 + (z021 + qnorm(1 - alpha/2)) / (1 - a21 * (z021 + qnorm(1 - alpha/2))))

  ##################for theta2
  alpha122 <- pnorm(z022 + (z022 + qnorm(alpha/2)) / (1 - a22 * (z022 + qnorm(alpha/2))))
  alpha222 <- pnorm(z022 + (z022 + qnorm(1 - alpha/2)) / (1 - a22 * (z022 + qnorm(1 - alpha/2))))

  ##################Calculate the lower and upper bound
  #### Compare to 0 and 1 (for p)!
  theta1_low <- theta1_bootstrap_ordered[as.integer(alpha11 * B)]
  theta1_up <- theta1_bootstrap_ordered[as.integer(alpha21 * B)]

  p_low <- p_bootstrap_ordered[as.integer(alpha1p * B)]
  p_low <- ifelse(p_low < 0, 0L, p_low)
  p_up <- p_bootstrap_ordered[as.integer(alpha2p * B)]
  p_up <- ifelse(p_up > 1, 1L, p_up)

  theta21_low <- theta21_bootstrap_ordered[as.integer(alpha121 * B)]
  theta21_up <- theta21_bootstrap_ordered[as.integer(alpha221 * B)]

  theta22_low <- theta22_bootstrap_ordered[as.integer(alpha122 * B)]
  theta22_up <- theta22_bootstrap_ordered[as.integer(alpha222 * B)]


  return(list(CI_theta_1 = c(theta1_low, theta1_up),
              CI_theta_21 = c(theta21_low, theta21_up),
              CI_theta_22 = c(theta22_low, theta22_up),
              CI_theta_p = c(p_low, p_up)))

}
