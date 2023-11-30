MLE_Geo<- function( data , n , tau , delta, theta21 , theta22 , p, maxit, tol, language, parallel, ncores ){

  # All Type I censoring




  ############################################################################
  ############# Part 1 first stress level
  ############################################################################

  # needs raw data

  # T1 <- sort(data)[sort(data)<tau1]
  # n1 <- length(T1)
  # T2 <- sort(data)[sort(data)>=tau1]

  ### q1 is the number of intervals under the first stress level
  q1 <- tau[1]/delta


  ### Count the number in each interval
  # n1j <- table(cut(T1, breaks = seq(0, tau1, delta)))
  # n1j <- table_factor_cpp(cut(T1, breaks = seq(0, tau1, delta)))

  n1j <- data[1:q1]
  n1 <- sum(n1j)

  if (n1 < 1){
    warning("No observation under the first stress level!")
    mle1 =NA
  }else{
    ### tau1j are the right bounds of intervals
    tau1j <- seq(delta, tau[1], delta)
    ### tau1j0 are the left bounds of intervals
    tau1j0 <- tau1j - delta

    ### Compute mle1
    ### Find the unique solution (it is proved to be unique) with uniroot function
    mle1 <- uniroot(func_theta1_int, c(1, 1000), n1j = n1j, n1 = n1, n=n, tau1=tau[1], tau1j=tau1j, tau1j0=tau1j0)$root
  }

  ############################################################################
  ############# Part 2 second stress level
  ############################################################################
  # n21 <- round(mean(p)*(n-n1))
  # n22 <- n-n1-n21
  # set.seed(1234)
  # T22 <- sample(T2, n22)
  # set.seed(1234)
  # T21 <-stl_sort(-theta21*log(1-runif(n21)) + tau1)
  # T2_combine <- stl_sort(c(T21, T22))
  # T2_combine
  # ### number of failures under the second stress level before the experiment stops
  #
  # # n2 <- sum(T2_combine < tau2) #not used in further calculations
  #
  # ### q2 is the number of intervals under the second stress level
  # q2 <- (tau2-tau1)/delta
  # j <- 1:q2
  # ### Count the number in each interval
  # # n2j
  # n2j <- table_factor_cpp(cut(T2_combine[T2_combine < tau2], breaks = seq(tau1, tau2, delta)))
  #
  # data <- rep(j-1, n2j)
  # data <- c(data, rep(q2, n-n1-n2))
  # # d <- as.numeric(data < q2) # input variable
  #
  # # supports both vectors and numeric values for p, theta21, theta22
  # parameter_starts <- expand.grid(p, theta21, theta22)
  # colnames(parameter_starts) <- c("omega1", "theta21", "theta22")

  ### q2 is the number of intervals under the second stress level
  q2 <- (tau[2]-tau[1])/delta

  if(q2 <= 3){
    warning("At least 3 different intervals in the second stress level are needed. The number of observations cannot perform an EM algorithm in a proper way.")
  }

  j <- 1:q2
  ### Count the number in each interval
  # if (n > sum(data)) {
  #   n2j <- data[-(1:q1)] # censored
  #   if(q2 > length(n2j)){
  #     n2j <- c(n2j, rep(0, q2-length(n2j)))
  #   }else if(q2 <- length(n2j)){
  #     n2j <- n2j[j]
  #   }
  # }else{n2j <- data[-c((1:q1),((q1+q2+1):length(data)))]}  #uncensored

  # Alternative
  if (n > sum(data)) { #data is censored with a tau2 that is smaller than the input tau2
    n2j <- data[-(1:q1)] # censored
    if(q2 > length(n2j)){n2j <- c(n2j, rep(0, q2-length(n2j)))}
  }else{n2j <- data[-c((1:q1),((q1+q2+1):length(data)))]}

  # # second alternative
  # n2j <- data[-c((1:q1),((q1+q2+1):length(data)))]
  # j <- 1:length(n2j)



  n2 <- sum(n2j)

  data_starts <- rep(j-1, n2j)
  data_starts <- c(data_starts, rep(q2, n-n1-n2))
  d <- as.numeric(data_starts < q2)

  parameter_starts <- data.frame(omega1=p, theta21=theta21, theta22)


  if(language == "CPP") {

    if (dim(parameter_starts)[1]==1){
      EM_mle <- EM_algorithm_interval_arma(ind =1, data = data_starts, N = maxit,
                                            delta = delta, d=d, parameter_starts=parameter_starts, q2=q2, tol = tol)
      # Remove index (unnecessary for one setup of parameters)
      # EM_mle$ind <-NULL
      EM_mle$results$ind <-NULL
      return(list(results=c(mle1=mle1, n1=n1, mle2=as.list(EM_mle), n2=n2, censored_rate=1-(n1+n2)/n), posterior=EM_mle$posterior))
    }else{

      if(parallel == TRUE){
      cl <- parallel::makeCluster(getOption("cl.cores", ncores))
      EM_mle <- parLapply(cl, 1:nrow(parameter_starts), EM_algorithm_interval_arma, data = data_starts, N = maxit,
                delta = delta, d=d, parameter_starts=parameter_starts, q2=q2, tol=tol)
      parallel::stopCluster(cl)

      }else {

      EM_mle <- lapply(1:nrow(parameter_starts), EM_algorithm_interval_arma, data = data_starts, N = maxit,
                          delta = delta, d=d, parameter_starts=parameter_starts, q2=q2, tol = tol)
      }

      Estimate_df <- data.frame(matrix(nrow = length(EM_mle), ncol = 8))
      colnames(Estimate_df) <- colnames(EM_mle[[1]]$results)
      posterior_l <- vector(mode = "list", length=length(EM_mle))
      for (i in 1:length(EM_mle)){
        # Estimate_df <- rbind.data.frame(Estimate_df, model_l7[[i]]$results)
        Estimate_df[i,] <- EM_mle[[i]]$results
        posterior_l[[i]] <- EM_mle[[i]]$posterior
      }

      # Estimate_df <- do.call("rbind.data.frame", EM_mle)


      mle_em <- find_max(Estimate_df)
      mle_em



      return(list(mle1=mle1, n1=n1, mle2=mle_em, n2=n2, censored_rate=1-(n1+n2)/n, estimates=Estimate_df, posterior=posterior_l))
    }

  }else{
    if (dim(parameter_starts)[1]==1){
      EM_mle <- EM_algorithm_interval(ind =1, data = data_starts, N = maxit,
                                           delta = delta, d=d, parameter_starts=parameter_starts, q2=q2, tol = tol)
      # Remove index (unnecessary for one setup of parameters)
      # EM_mle$ind <-NULL
      EM_mle$results$ind <-NULL
      return(list(results=c(mle1=mle1, n1=n1, mle2=as.list(EM_mle), n2=n2, censored_rate=1-(n1+n2)/n), posterior=EM_mle$posterior))
    }else{

      if(parallel == TRUE){
        cl <- parallel::makeCluster(getOption("cl.cores", ncores))
        EM_mle <- parLapply(cl, 1:nrow(parameter_starts), EM_algorithm_interval, data = data_starts, N = maxit,
                            delta = delta, d=d, parameter_starts=parameter_starts, q2=q2, tol = tol)
        parallel::stopCluster(cl)

      }else {

        EM_mle <- lapply(1:nrow(parameter_starts), EM_algorithm_interval, data = data_starts, N = maxit,
                         delta = delta, d=d, parameter_starts=parameter_starts, q2=q2, tol = tol)
      }

      # Estimate_df <- do.call("rbind.data.frame", EM_mle)

      Estimate_df <- data.frame(matrix(nrow = length(EM_mle), ncol = 8))
      colnames(Estimate_df) <- colnames(EM_mle[[1]]$results)
      posterior_l <- vector(mode = "list", length=length(EM_mle))
      for (i in 1:length(EM_mle)){
        # Estimate_df <- rbind.data.frame(Estimate_df, model_l7[[i]]$results)
        Estimate_df[i,] <- EM_mle[[i]]$results
        posterior_l[[i]] <- EM_mle[[i]]$posterior
      }


      mle_em <- find_max(Estimate_df)
      mle_em



      return(list(mle1=mle1, n1=n1, mle2=mle_em, n2=n2, censored_rate=1-(n1+n2)/n, estimates=Estimate_df, posterior=posterior_l))
    }
  }
}
