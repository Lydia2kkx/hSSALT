
bootstrap_continuous <-function(data, n, parameter_starts, censoring, tau , r, B,
                       p, maxit, tol, language, parallel, ncores){

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
    n21 <- rbinom(1, (n-n1), p)
    T21 <- stl_sort(-theta21*log(1-runif(n21)) + tau[1])
    n22 <- n-n1-n21
    T22 <- stl_sort(-theta22*log(1-runif(n22)) + tau[1])
    T2 <- stl_sort(c(T21, T22))
    n2 <- sum(T2 < tau[2])
    d <- c(rep(1, n2), rep(0, n-n1-n2))
    t22 <- apply(cbind(T2, tau[2]), 1, min) # observed or censored data



    # model_list <- lapply(1:nrow(parameter_starts), EM_algorithm_censored, data = t22 - tau[1])


    # parameter_starts <- data.frame(omega1=p, theta21=theta21, theta22)


    if(language == "CPP") {

      if (dim(parameter_starts)[1]==1){
        EM_mle <- EM_algorithm_censored_arma(ind =1, data=t22-tau[1], d=d, N=maxit, parameter_starts = parameter_starts, tol = tol)
        Estimate_df <- EM_mle$results

      }else if(parallel == TRUE){
        cl <- parallel::makeCluster(getOption("cl.cores", ncores))
        EM_mle <- parLapply(cl, 1:nrow(parameter_starts), EM_algorithm_censored_arma, data = t22 - tau[1], d=d, N=maxit,
                            parameter_starts = parameter_starts, tol = tol)
        parallel::stopCluster(cl)
        # Estimate_df <- do.call("rbind.data.frame", EM_mle)
        Estimate_df <- data.frame(matrix(nrow = length(EM_mle), ncol = 8))
        colnames(Estimate_df) <- colnames(EM_mle[[1]]$results)
        posterior_l <- vector(mode = "list", length=length(EM_mle))
        for (i in 1:length(EM_mle)){
          # Estimate_df <- rbind.data.frame(Estimate_df, model_l7[[i]]$results)
          Estimate_df[i,] <- EM_mle[[i]]$results
          posterior_l[[i]] <- EM_mle[[i]]$posterior
        }

      }else{
        EM_mle <- lapply(1:nrow(parameter_starts), EM_algorithm_censored_arma, data = t22 - tau[1], d=d, N=maxit,
                         parameter_starts = parameter_starts, tol = tol)
        # Estimate_df <- do.call("rbind.data.frame", EM_mle)
        Estimate_df <- data.frame(matrix(nrow = length(EM_mle), ncol = 8))
        colnames(Estimate_df) <- colnames(EM_mle[[1]]$results)
        posterior_l <- vector(mode = "list", length=length(EM_mle))
        for (i in 1:length(EM_mle)){
          # Estimate_df <- rbind.data.frame(Estimate_df, model_l7[[i]]$results)
          Estimate_df[i,] <- EM_mle[[i]]$results
          posterior_l[[i]] <- EM_mle[[i]]$posterior
        }
      }


    }else{

      if (dim(parameter_starts)[1]==1){
        EM_mle <- EM_algorithm_censored(ind =1, data=t22-tau[1], d=d, N=maxit, parameter_starts = parameter_starts, tol)
        Estimate_df <- EM_mle$results

      }else if(parallel == TRUE){
        cl <- parallel::makeCluster(getOption("cl.cores", ncores))
        EM_mle <- parLapply(cl, 1:nrow(parameter_starts), EM_algorithm_censored, data = t22 - tau[1], d=d, N=maxit,
                            parameter_starts = parameter_starts, tol)
        parallel::stopCluster(cl)
        Estimate_df <- data.frame(matrix(nrow = length(EM_mle), ncol = 8))
        colnames(Estimate_df) <- colnames(EM_mle[[1]]$results)
        posterior_l <- vector(mode = "list", length=length(EM_mle))
        for (i in 1:length(EM_mle)){
          # Estimate_df <- rbind.data.frame(Estimate_df, model_l7[[i]]$results)
          Estimate_df[i,] <- EM_mle[[i]]$results
          posterior_l[[i]] <- EM_mle[[i]]$posterior
        }

      }else{
        EM_mle <- lapply(1:nrow(parameter_starts), EM_algorithm_censored, data = t22 - tau[1], d=d, N=maxit,
                         parameter_starts = parameter_starts, tol)
        Estimate_df <- data.frame(matrix(nrow = length(EM_mle), ncol = 8))
        colnames(Estimate_df) <- colnames(EM_mle[[1]]$results)
        posterior_l <- vector(mode = "list", length=length(EM_mle))
        for (i in 1:length(EM_mle)){
          # Estimate_df <- rbind.data.frame(Estimate_df, model_l7[[i]]$results)
          Estimate_df[i,] <- EM_mle[[i]]$results
          posterior_l[[i]] <- EM_mle[[i]]$posterior
        }
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
# bootstrap_distri <- cbind.data.frame(i = 1:nrow(max_loglik_df), max_loglik_df)

return(cbind.data.frame(i = 1:nrow(max_loglik_df), max_loglik_df))

}
