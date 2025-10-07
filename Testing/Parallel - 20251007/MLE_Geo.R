MLE_Geo<- function(data, n, tau, delta, theta21, theta22, p, maxit, tol, language, parallel, ncores){

  # All Type I censoring

  ############################################################################
  ############# Part 1 first stress level
  ############################################################################

  ### q1 is the number of intervals under the first stress level
  q1 <- tau[1]/delta


  ### Count the number in each interval
  n1j <- data[1:q1]
  n1 <- sum(n1j)

  if (n1 < 1){
    warning("No observation under the first stress level!")
    mle1 = NA
  }else{
    ### tau1j are the right bounds of intervals
    tau1j <- seq(delta, tau[1], delta)
    ### tau1j0 are the left bounds of intervals
    tau1j0 <- tau1j - delta

    ### Compute mle1
    ### Find the unique solution (it is proved to be unique) with uniroot function
    mle1 <- uniroot(func_theta1_int, c(1, 1000), n1j = n1j, n1 = n1, n = n, 
                    tau1 = tau[1], tau1j=tau1j, tau1j0=tau1j0)$root
  }

  ### q2 is the number of intervals under the second stress level
  q2 <- (tau[2]-tau[1])/delta

  if(q2 <= 3){
    warning("At least 3 different intervals in the second stress level are needed. The number of observations cannot perform an EM algorithm in a proper way.")
  }

  j <- 1:q2

  n2j <- data[(q1+1):(q1+q2)] #Avner: Now works even if the inputted data is the full data (where the previous version failed)
  if(q1+q2 > length(data)){n2j[(length(data)-q1+1):q2] <- 0} #Avner: If the inputted data does not cover all of n2 then fill with zeros
  n2 <- sum(n2j)
  
  data_starts <- rep(j-1, n2j)
  data_starts <- c(data_starts, rep(q2, n-n1-n2))
  d <- as.numeric(data_starts < q2)

  parameter_starts <- data.frame(omega1=p, theta21=theta21, theta22=theta22)

  if(language == "CPP") {

    if (dim(parameter_starts)[1]==1){
      EM_mle <- EM_algorithm_interval_arma(ind =1, data = data_starts, N = maxit,
                                            delta = delta, d=d, parameter_starts=parameter_starts, q2=q2, tol = tol)
    }else{

      if(parallel == TRUE){
      cl <- parallel::makeCluster(getOption("cl.cores", ncores))
      parallel::clusterExport(cl, varlist = c("mysum_cpp", "EM_algorithm_interval_arma"))
      EM_mle <- parallel::parLapply(cl, 1:nrow(parameter_starts), EM_algorithm_interval_arma, data = data_starts, N = maxit,
                delta = delta, d=d, parameter_starts=parameter_starts, q2=q2, tol=tol)
      parallel::stopCluster(cl)

      }else {

      EM_mle <- lapply(1:nrow(parameter_starts), EM_algorithm_interval_arma, data = data_starts, N = maxit,
                          delta = delta, d=d, parameter_starts=parameter_starts, q2=q2, tol = tol)
      }
    }

  }else{
    if (dim(parameter_starts)[1]==1){

      EM_mle <- EM_algorithm_interval(ind =1, data = data_starts, N = maxit,
                                           delta = delta, d=d, parameter_starts=parameter_starts, q2=q2, tol = tol)
    }else{

      if(parallel == TRUE){
        cl <- parallel::makeCluster(getOption("cl.cores", ncores))
        parallel::clusterExport(cl, varlist = c("sum_finite", "EM_algorithm_interval"))
        EM_mle <- parallel::parLapply(cl, 1:nrow(parameter_starts), EM_algorithm_interval, data = data_starts, N = maxit,
                            delta = delta, d=d, parameter_starts=parameter_starts, q2=q2, tol = tol)
        parallel::stopCluster(cl)

      }else {

        EM_mle <- lapply(1:nrow(parameter_starts), EM_algorithm_interval, data = data_starts, N = maxit,
                         delta = delta, d=d, parameter_starts=parameter_starts, q2=q2, tol = tol)
      }
    }
  }
  
  if (dim(parameter_starts)[1]==1){
    # Remove index (unnecessary for one setup of parameters)
    # EM_mle$ind <-NULL
    EM_mle$results$ind <-NULL
    
    if (!is.na(EM_mle$results$theta22) && !is.na(mle1) && mle1 < EM_mle$results$theta22) {
      mle1 = NA
      EM_mle$results$theta22 = NA
      warning("Mean lifetime in first stress level must be greater than mean lifetimes in the second stress level.") 
    }
    
    #Yao: added check theta21 = 0.
    #Yao: I think all these checks can be removed in the end to check only once. This change will be done later.
    if (EM_mle$results$theta21 == 0) {
      warning("The group with the smaller mean lifetime at the second stress level is estimated to be 0") 
    }
    
    #Check homogeneity
    if (!is.na(EM_mle$results$theta22)) {
      info <- ifelse(EM_mle$results$theta22/EM_mle$results$theta21 < 1.05, "homogeneous", "heterogeneous")
    } else {
      info <- NA
    }
    
    mle <- cbind(EM_mle$results[,1:2],theta1=mle1,EM_mle$results[,3:4])
    output <- list(n1=n1, n2=n2, mle=mle, loglik=EM_mle$results$loglik, iteration=EM_mle$results$iteration, message=EM_mle$results$message, info=info,censored_rate=1-(n1+n2)/n, posterior=EM_mle$posterior)
  } else {
    Estimate_df <- data.frame(matrix(nrow = length(EM_mle), ncol = 8))
    colnames(Estimate_df) <- colnames(EM_mle[[1]]$results)
    posterior_l <- vector(mode = "list", length=length(EM_mle))
    for (i in 1:length(EM_mle)){
      Estimate_df[i,] <- EM_mle[[i]]$results
      posterior_l[[i]] <- EM_mle[[i]]$posterior
    }
    
    mle_em <- find_max(Estimate_df)
    # Check whether EM algorithm converged
    if (nrow(mle_em)>0 && !is.na(mle1)) {
      posterior_l <- EM_mle[[as.numeric(rownames(Estimate_df[which.max(Estimate_df$loglik),]))]]$posterior
      # Check whether mle1 > theta 22
      if (mle1 < mle_em$theta22) {
        mle1 = NA
        mle_em$theta22 = NA
        warning("Mean lifetime in first stress level must be greater than mean lifetimes in the second stress level.")  
      }
    } else {
      posterior_l <- NA
    }
    mle <- cbind(mle_em[,1:2],theta1=mle1,mle_em[,3:4])
    
    #Check homogeneity
    if (!is.na(mle$theta22)) {
      info <- ifelse(mle$theta22/mle$theta21 < 1.05, "homogeneous", "heterogeneous")
    } else {
      info <- NA
    }
    
    output <- list(n1=n1, n2=n2, mle=mle,loglik=mle_em$loglik,iteration=mle_em$iteration,message=mle_em$message, info=info, censored_rate=1-(n1+n2)/n, posterior=posterior_l)
  }
    
  class(output) <- "hSSALTMLE"
  return(output)
}
