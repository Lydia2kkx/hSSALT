
CIbca_hSSALT<- function(data , n , censoring , tau , r, monitoring , alpha , B,
                        theta1, theta21 , theta22 , p , maxit , tol , language , parallel , ncores ){


  ########## Continuous ##########

  # T1 <- data[data<tau[1]]
  # n1 <- length(T1)
  # T2 <- data[data>=tau[1]]
  # n2 <- length(T2)
  # values of real data
  same_threshold <- 0.05
  p_threshold <- 0.01
  p_adjust <- TRUE ##adjust p with same results to 0.5
  data_adjust <- TRUE ##discard same results of theta21 and theta22
  mle1<-theta1

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
          Estimate_df <- EM_mle

        }else if(parallel == TRUE){
          cl <- parallel::makeCluster(getOption("cl.cores", ncores))
          EM_mle <- parLapply(cl, 1:nrow(parameter_starts), EM_algorithm_censored_arma, data = t22 - tau[1], d=d, N=maxit,
                              parameter_starts = parameter_starts, tol = tol)
          parallel::stopCluster(cl)
          Estimate_df <- do.call("rbind.data.frame", EM_mle)

        }else{
          EM_mle <- lapply(1:nrow(parameter_starts), EM_algorithm_censored_arma, data = t22 - tau[1], d=d, N=maxit,
                           parameter_starts = parameter_starts, tol = tol)
          Estimate_df <- do.call("rbind.data.frame", EM_mle)
        }


      }else{

        if (dim(parameter_starts)[1]==1){
          EM_mle <- EM_algorithm_censored(ind =1, data=t22-tau[1], d=d, N=maxit, parameter_starts = parameter_starts, tol)
          Estimate_df <- EM_mle

        }else if(parallel == TRUE){
          cl <- parallel::makeCluster(getOption("cl.cores", ncores))
          EM_mle <- parLapply(cl, 1:nrow(parameter_starts), EM_algorithm_censored, data = t22 - tau[1], d=d, N=maxit,
                              parameter_starts = parameter_starts, tol)
          parallel::stopCluster(cl)
          Estimate_df <- do.call("rbind.data.frame", EM_mle)

        }else{
          EM_mle <- lapply(1:nrow(parameter_starts), EM_algorithm_censored, data = t22 - tau[1], d=d, N=maxit,
                           parameter_starts = parameter_starts, tol)
          Estimate_df <- do.call("rbind.data.frame", EM_mle)
        }

      }


      Estimate_df$theta1 <- theta1_hat
      Estimate_df$n2 <- n2
      Estimate_df$p1 <- ifelse(Estimate_df$theta22 >= theta1_hat, NA, Estimate_df$p1)
      Estimate_df$info <- ifelse((Estimate_df$p1 < p_threshold) | (Estimate_df$p2 < p_threshold) | (Estimate_df$theta22 - Estimate_df$theta21 <= same_threshold*Estimate_df$theta22),
                                 "homo", ifelse(is.na(Estimate_df$p1), NA, "heter"))
      Estimate_df$mle2_classi <- (sum(T2[1:n2] - tau[1]) + (n - n1 - n2) * (tau[2] - tau[1]))/ n2
      if(p_adjust){
        Estimate_df$p1[Estimate_df$info == "homo"] <- Estimate_df$p2[Estimate_df$info == "homo"] <- 0.5
      }
      if(data_adjust){
        Estimate_df$p1[Estimate_df$info == "homo"] <- NA
      }
      if(all(is.na(Estimate_df$p1))){
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

  ############Boostrap BCa CIs
  ### Compute the bias correction
  ### order the bootstrap estimates
  theta1_bootstrap_ordered <- stl_sort(bootstrap_distri$theta1)
  p_bootstrap_ordered <- stl_sort(bootstrap_distri$p1)
  theta21_bootstrap_ordered <- stl_sort(bootstrap_distri$theta21)
  theta22_bootstrap_ordered <- stl_sort(bootstrap_distri$theta22)

  z01 <- qnorm(sum(theta1_bootstrap_ordered < mle1) / B)
  z0p <- qnorm(sum(p_bootstrap_ordered < p) / B)
  z021 <- qnorm(sum(theta21_bootstrap_ordered < theta21) / B)
  z022 <- qnorm(sum(theta22_bootstrap_ordered < theta22) / B)

  ### BCa method is computed based on the sample data
  ### Compute the acceleration ai
  ####Calculation of a1
  ### T1 is the failures under the first stress level (simulated or real ??)
  ### Again, take Type-I as an example
  mle1ij <- c()
  for (i in 1 : n1) {
    T1_new <- T1[-i]
    mle_1new <- (sum(T1_new) + (n - n1) * tau[1]) / (n1 - 1)
    mle1ij <- c(mle1ij, mle_1new)
  }

  mle1ij <- c(mle1ij, rep((sum(T1) + (n - 1 - n1) * tau[1]) / n1, n2))
  mle1i <- mean(mle1ij)
  a1 <- sum((mle1i - mle1ij)^3) / (6 * (sum((mle1i - mle1ij)^2))^1.5)

  ####Calculation of a2s
  ### T2 is the failures under the second stress level
  Re <- list()
  for (i in 1 : (length(T2) - (n-n1-n2))) {
    T2_new <- T2[-i]
    d <- 1*(T2_new <= tau[2])
    t22_new <- apply(cbind(T2_new, tau[2]), 1, min) # observed or censored data
    model_list_new <- lapply(1:nrow(parameter_starts), EM_algorithm_censored, data = t22_new - tau[1])
    Estimate_df_new <- do.call("rbind.data.frame", model_list_new)
    Re[[i]] <- Estimate_df_new
  }

  ### Again here, if the initial values are vectors, use find_max
  ### Otherwise, simply skip
  stab_re <- lapply(Re, check_stability)
  stab_re_df <- do.call("rbind.data.frame", stab_re)
  max_loglik_re <- lapply(Re, find_max)
  max_loglik_re_df <- do.call("rbind.data.frame", max_loglik_re)


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
  p_up <- p_bootstrap_ordered[as.integer(alpha2p * B)]


  theta21_low <- theta21_bootstrap_ordered[as.integer(alpha121 * B)]
  theta21_up <- theta21_bootstrap_ordered[as.integer(alpha221 * B)]


  theta22_low <- theta22_bootstrap_ordered[as.integer(alpha122 * B)]
  theta22_up <- theta22_bootstrap_ordered[as.integer(alpha222 * B)]


  return(list(CI_theta_1 = c(theta1_low, theta1_up),
              CI_theta_21 = c(theta21_low, theta21_up),
              CI_theta_22 = c(theta22_low, theta22_up),
              CI_theta_p = c(p_low, p_up)
  ))


}
