MLE_Exp <- function( data , n , censoring , tau , r , theta21 , theta22 , p, maxit, tol, language, parallel, ncores ){




  # n <- length(data)

  # in type 1 censoring tau=tau1 and r=tau2

  ############################################################################
  ############# Part 1 first stress level
  ############################################################################

  # T1 <- sort(data)[sort(data)<tau]
  # n1 <- length(T1)
  # T2 <- sort(data)[sort(data)>=tau]

  T1 <- data[data<tau[1]]
  n1 <- length(T1)
  T2 <- data[data>=tau[1]]

  if (n1 < 1){
    warning("No observation under the first stress level!")
    mle1 =NA
  }else{
  ### Find the unique solution (it is proved to be unique) with uniroot function
  mle1 <- (sum(T1) + tau[1]*(n-n1))/n1
  }

  ############################################################################
  ############# Part 2 second stress level
  ############################################################################


  n2 <- length(T2)
  if (n2 <=2){
    warning("At least 3 observations in the second stress level are needed. The number of observations cannot perform an EM algorithm in a proper way.")
  }

  # n21 <- round(mean(p)*(n-n1))
  # n22 <- n-n1-n21
  # set.seed(1234)
  # T22 <- sample(T2, n22)
  # set.seed(1234)
  # T21 <- stl_sort(-theta21*log(1-runif(n21)) + tau)
  # T2_combine <- stl_sort(c(T21, T22))

  #
  # if (censoring ==1){
  #
  #   ### number of failures under the second stress level before the experiment stops
  #   # n2 <- sum(T2_combine < tau2) # not relevant
  #   # cs <- T2_combine[r-n1]
  #   # d <- 1*(T2_combine < r)              # r = tau2
  #   # t22 <- T2_combine # observed or censored data
  #   # t22[which(t22>=r)] = r #set data after censoring tau2 (r) to tau 2
  #
  #   d <- 1*(T2 < tau[2])              # r = tau2
  #   t22 <- T2 # observed or censored data
  #   t22[which(t22>=tau[2])] = tau[2] #set data after censoring tau2 (r) to tau 2
  #
  # }
  #
  # if (censoring ==2){
  #   # cs <- T2_combine[r-n1]
  #   # d <- 1*(T2_combine <= cs)              # censoring indicator
  #   # t22 <- apply(cbind(T2_combine, cs), 1, min) # observed or censored data
  #
  #   cs <- T2[r-n1]
  #   d <- 1*(T2 <= cs)              # censoring indicator
  #   t22 <- apply(cbind(T2, cs), 1, min) # observed or censored data
  #
  # }


  # retrieves same result for censored and uncensored data

  n_c <- length(data)

  if (censoring ==1){

    if (n>n_c){  #censored data
      t22 <- c(T2,rep(tau[2], n-n_c))
      d <- 1*(t22 < tau[2])
    }else{       #uncensored data
      d <- 1*(T2 < tau[2])              # r = tau2
      t22 <- T2 # observed or censored data
      t22[which(t22>=tau[2])] = tau[2] #set data after censoring tau2 to tau 2
    }
  }

  if (censoring ==2){
    cs <- T2[r-n1]

    if (n>n_c){
      t22 <- c(T2, rep(cs,n-n_c))
      d <- c(rep(1,n_c-n1), rep(0,n-n_c))
    }else{
      d <- 1*(T2 <= cs)              # censoring indicator
      t22 <- apply(cbind(T2, cs), 1, min) # observed or censored data
    }
  }



  # supports both vectors and numeric values for p, theta21, theta22
  # parameter_starts <- expand.grid(p, theta21, theta22)
  # colnames(parameter_starts) <- c("omega1", "theta21", "theta22")

  parameter_starts <- data.frame(omega1=p, theta21=theta21, theta22)


  if(language == "CPP") {

    if (dim(parameter_starts)[1]==1){
      EM_mle <- EM_algorithm_censored_arma(ind =1, data=t22-tau[1], d=d, N=maxit, parameter_starts = parameter_starts, tol = tol)
      # Remove index (unnecessary for one setup of parameters)
      # EM_mle$ind <-NULL
      EM_mle$results$ind <-NULL
      # return(c(mle1=mle1, n1=n1, mle2=as.list(EM_mle), n2=n2, censored_rate=1-n_c/n ))
      return(list(results=c(mle1=mle1, n1=n1, mle2=as.list(EM_mle$results), n2=n2, censored_rate=1-n_c/n ), posterior=EM_mle$posterior))
    }else{

      if(parallel == TRUE){
        cl <- parallel::makeCluster(getOption("cl.cores", ncores))
        EM_mle <- parLapply(cl, 1:nrow(parameter_starts), EM_algorithm_censored_arma, data = t22 - tau[1], d=d, N=maxit,
                            parameter_starts = parameter_starts, tol = tol)
        parallel::stopCluster(cl)

      }else {

        EM_mle <- lapply(1:nrow(parameter_starts), EM_algorithm_censored_arma, data = t22 - tau[1], d=d, N=maxit,
                         parameter_starts = parameter_starts, tol = tol)
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



      return(list(mle1=mle1, n1=n1, mle2=mle_em, n2=n2, censored_rate=1-n_c/n, estimates=Estimate_df, posterior=posterior_l))
    }

  }else{

    if (dim(parameter_starts)[1]==1){
      EM_mle <- EM_algorithm_censored(ind =1, data=t22-tau[1], d=d, N=maxit, parameter_starts = parameter_starts, tol)
      # Remove index (unnecessary for one setup of parameters)
      # EM_mle$ind <-NULL
      EM_mle$results$ind <-NULL
      # return(c(mle1=mle1, n1=n1, mle2=as.list(EM_mle), n2=n2, censored_rate=1-n_c/n ))
      return(list(results=c(mle1=mle1, n1=n1, mle2=as.list(EM_mle$results), n2=n2, censored_rate=1-n_c/n ), posterior=EM_mle$posterior))
    }else{

      if(parallel == TRUE){
        cl <- parallel::makeCluster(getOption("cl.cores", ncores))
        EM_mle <- parLapply(cl, 1:nrow(parameter_starts), EM_algorithm_censored, data = t22 - tau[1], d=d, N=maxit,
                            parameter_starts = parameter_starts, tol)
        parallel::stopCluster(cl)

      }else {

        EM_mle <- lapply(1:nrow(parameter_starts), EM_algorithm_censored, data = t22 - tau[1], d=d, N=maxit,
                         parameter_starts = parameter_starts, tol)
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



      return(list(mle1=mle1, n1=n1, mle2=mle_em, n2=n2, censored_rate=1-n_c/n, estimates=Estimate_df, posterior=posterior_l))
    }



  }

}

