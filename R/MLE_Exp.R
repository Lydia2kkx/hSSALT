MLE_Exp <- function(data, n, censoring, tau, r, theta21, theta22, p, maxit, tol, language){
  
  ############################################################################
  ############# Part 1 first stress level
  ############################################################################
  data <- sort(data)
  T1 <- data[data < tau[1]]
  n1 <- length(T1)
  T2 <- data[data >= tau[1]]
  
  if (n1 < 1){
    warning("No observation under the first stress level!")
    mle1 <- NA
  }else{
    mle1 <- (sum(T1) + tau[1]*(n-n1))/n1
  }
  
  ############################################################################
  ############# Part 2 second stress level
  ############################################################################
  
  n_c <- length(data)
  
  if (censoring == 1){
    n2 <- length(T2[T2 <= tau[2]])
    
    if (n > n_c){
      t22 <- c(T2, rep(tau[2], n-n_c))
      d <- 1*(t22 < tau[2])
    }else{
      d <- 1*(T2 < tau[2])
      t22 <- T2
      t22[which(t22 >= tau[2])] = tau[2]
    }
  }
  
  if (censoring == 2){
    n2 <- r-n1
    cs <- T2[r-n1]
    
    if (n>n_c){
      t22 <- c(T2, rep(cs, n-n_c))
      d <- c(rep(1, n_c-n1), rep(0, n-n_c))
    }else{
      d <- 1*(T2 <= cs)
      t22 <- apply(cbind(T2, cs), 1, min)
    }
  }
  
  if (n2 <= 2){
    warning("At least 3 observations in the second stress level are needed. The number of observations cannot perform an EM algorithm in a proper way.")
  }
 
  parameter_starts <- data.frame(omega1 = p, theta21 = theta21, theta22 = theta22)

  if(language == "CPP") {
    
    if (dim(parameter_starts)[1]==1){
      EM_mle <- EM_algorithm_censored_arma(ind =1, data=t22-tau[1], d=d, N=maxit, parameter_starts = parameter_starts, tol = tol)
    }else{
      EM_mle <- lapply(1:nrow(parameter_starts), EM_algorithm_censored_arma, data = t22 - tau[1], d=d, N=maxit,
                         parameter_starts = parameter_starts, tol = tol)
    }
    
  }else{
    
    if (dim(parameter_starts)[1] == 1){
      EM_mle <- EM_algorithm_censored(ind = 1, data = t22-tau[1], d = d, N = maxit, 
                                      parameter_starts = parameter_starts, tol = tol)
    }else{
      EM_mle <- lapply(1:nrow(parameter_starts), EM_algorithm_censored, data = t22 - tau[1], d=d, N=maxit,
                         parameter_starts = parameter_starts, tol)
    }
  }
  
  if (dim(parameter_starts)[1]==1){
    EM_mle$results$ind <- NULL
    
    ###Check whether has converged and mle1 > theta22
    if (!is.na(EM_mle$results$theta22) && !is.na(mle1) && mle1 < EM_mle$results$theta22) {
      mle1 = NA
      EM_mle$results$theta22 = NA
      warning("Mean lifetime in first stress level must be greater than mean lifetimes in the second stress level.") 
    }
    
    ###Check homogeneity
    if (!is.na(EM_mle$results$theta22)) {
      info <- ifelse(EM_mle$results$theta22/EM_mle$results$theta21 < 1.05, "homogeneous", "heterogeneous")
    } else {
      info <- NA
    }
    
    mle <- cbind(EM_mle$results[, 1:2], theta1 = mle1, EM_mle$results[, 3:4])
    output <- list(n1 = n1, n2 = n2, mle = mle, loglik = EM_mle$results$loglik, 
                   iteration = EM_mle$results$iteration, message = EM_mle$results$message, 
                   info = info, censored_rate = 1-(n1+n2)/n, posterior = EM_mle$posterior)
  } else {
    Estimate_df <- data.frame(matrix(nrow = length(EM_mle), ncol = 7))
    colnames(Estimate_df) <- colnames(EM_mle[[1]]$results)
    for (i in 1:length(EM_mle)){
      Estimate_df[i,] <- EM_mle[[i]]$results
    }
    
    mle_em <- find_max(Estimate_df)
    
    ###Check whether EM algorithm converged
    if (nrow(mle_em)>0 && !is.na(mle1)) {
      posterior_l <- EM_mle[[as.numeric(rownames(Estimate_df[which.max(Estimate_df$loglik),]))]]$posterior
      
      ###Check whether mle1 > theta 22
      if (mle1 < mle_em$theta22) {
        mle1 = NA
        mle_em$theta22 = NA
        warning("Mean lifetime in first stress level must be greater than mean lifetimes in the second stress level.")  
      }
      
    } else {
      posterior_l <- NA
    }
    
    ###Check homogeneity
    if (!is.na(mle_em$theta22)) {
      if (mle_em$theta22/mle_em$theta21 < 1.05) {
        warning("The dataset appears homogeneous!")
      }
    }
    
    mle <- cbind(mle_em[, 1:2], theta1 = mle1, mle_em[, 3:4])
    
    ###Check homogeneity
    if (!is.na(mle$theta22)) {
      info <- ifelse(mle$theta22/mle$theta21 < 1.05, "homogeneous", "heterogeneous")
    } else {
      info <- NA
    }
    
    output <- list(n1 = n1, n2 = n2, mle = mle, loglik = mle_em$loglik, 
                   iteration = mle_em$iteration, message = mle_em$message, 
                   info = info, censored_rate = 1-(n1+n2)/n, posterior = posterior_l)
  }
  
  class(output) <- "hSSALTMLE"
  return(output)
}