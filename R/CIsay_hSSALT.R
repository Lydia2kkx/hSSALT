

CIsay_hSSALT <- function(data, n, censoring, tau , r, monitoring, delta, alpha, theta1, p, theta21, theta22){
  #Yao: please also consider the data as full or censored data, adjust the code accordingly. 
  # Interval case
  if(monitoring=="interval"){
    
    # get interval bounds from data
    
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
    
    # if (n > sum(data)) { #data is censored with a tau2 that is smaller than the input tau2
    #   n2j <- data[-(1:q1)] # censored
    #   if(q2 > length(n2j)){n2j <- c(n2j, rep(0, q2-length(n2j)))}
    # }else{n2j <- data[-c((1:q1),((q1+q2+1):length(data)))]}
    
    #Avner: My version of the above code - Added after testing 20250901
    n2j <- data[(q1+1):(q1+q2)] 
    if(q1+q2 > length(data)){n2j[(length(data)-q1+1):q2] <- 0}
    n2 <- sum(n2j)
    
    Secondderivative_O11 <- function(theta1, tau1j, tau1j0, n1j, n1, tau1){
      -2/(theta1^3)*sum(n1j*(tau1j0*exp(-tau1j0/theta1)-tau1j*exp(-tau1j/theta1))/(exp(-tau1j0/theta1)-exp(-tau1j/theta1))) + 
        1/(theta1^2)*sum(n1j*((tau1j0/theta1)^2*exp(-tau1j0/theta1)-(tau1j/theta1)^2*exp(-tau1j/theta1))/(exp(-tau1j0/theta1)-exp(-tau1j/theta1)))-
        1/(theta1^2)*sum(n1j*((tau1j0/theta1*exp(-tau1j0/theta1)-tau1j/theta1*exp(-tau1j/theta1))/(exp(-tau1j0/theta1)-exp(-tau1j/theta1)))^2)-
        2*tau1*(n-n1)/(theta1^3)
    }
    O11 <- -Secondderivative_O11(theta1 = theta1, tau1j, tau1j0, n1j, n1, tau1 = tau[1])
    V_theta1 <- 1/O11
    
    # alpha is the significance level 
    theta1_approxCI_low <- theta1 - qnorm(1-alpha/2)*sqrt(V_theta1)
    theta1_approxCI_low <- ifelse(theta1_approxCI_low < 0, 0L, theta1_approxCI_low)
    theta1_approxCI_up <- theta1 + qnorm(1-alpha/2)*sqrt(V_theta1)
    
    nf <- n1+n2
    
    loglik_i_interval <- function(x){
      theta21 <- x[1]
      theta22 <- x[2]
      p <- x[3]
      
      value1 <- p*(exp(-(tau2j0-tau[1])/theta21)-exp(-(tau2j-tau[1])/theta21))
      value2 <- (1-p)*(exp(-(tau2j0-tau[1])/theta22)-exp(-(tau2j-tau[1])/theta22))
      
      sum(n2j*log(value1+value2)) +
        (n-nf)*log(p*exp(-(tau[2]-tau[1])/theta21)+(1-p)*exp(-(tau[2]-tau[1])/theta22))
      
    }
    
    x <- c(theta21,theta22,p)
    hess <- numDeriv::hessian(func = loglik_i_interval, x = x)
  
    
  }else{############ continuous ###############
    
    mle1 <- theta1
    
    if(censoring==1){
      
      # Theta1, Type I
      
      #Avner: Changed from 50 to 30 to check stability
      if(n>25){
        suppressWarnings({
        ###The case with a relatively large sample size (e.g., n = 50)
        ###High precision computation is needed
        T1_simu <- Rmpfr::mpfr(1, 256)*data
        n1 <- sum(T1_simu <= tau[1])
        T1 <- sort(T1_simu[T1_simu <= tau[1]])
        mle1 <- Rmpfr::mpfr(sum(T1, (n-n1)*tau[1]) / n1, 256)
        
        p1 <- 1 - exp(-tau[1]/mle1)
        p2 <- exp(-tau[1]/mle1) * (1 - p1*exp(-(tau[2]-tau[1])/theta21) - (1-p1)*exp(-(tau[2]-tau[1])/theta22))
        p3 <- 1 - p1 - p2
        Cn <- 1/(1 - (1-p1)^n - (1-p2)^n + p3^n)
        
        # add mle1, tau
        bias_part_Mpfr <- function(n){
          y <- 0
          
          for (i in 1:(n-1)) {
            for (k in 0:i) {
              n <- Rmpfr::mpfr(n, 256)
              C_ik <- (-1)^k * Rmpfr::chooseMpfr(n, i) * Rmpfr::chooseMpfr(i, k) * ((1 - p1)^(n-i) - p3^(n-i)) * (1-p1)^k
              tau_ik <- tau[1] / i * (n - i + k)
              re_ik <- C_ik *tau_ik
              y <- y + re_ik
            }
          }
          return(y)
        }
        
        V1 <- mle1^2/n1
        # V1
        
        ####################################################################################
        ######Censored
        bias <- bias_part_Mpfr(n)
        
        theta1_approxCI_low <- mle1 - Cn * bias - qnorm(1-alpha/2)*sqrt(V1)
        theta1_approxCI_low <- ifelse(theta1_approxCI_low < 0, 0L, theta1_approxCI_low)
        #Avner: Reverted to previous version because new version caused errors.
        theta1_approxCI_low <- sapply(theta1_approxCI_low[1], Rmpfr::asNumeric)
        #theta1_approxCI_low <- Rmpfr::asNumeric(theta1_approxCI_low)
        theta1_approxCI_up <- mle1 - Cn * bias + qnorm(1-alpha/2)*sqrt(V1)
        #Avner: Same as above
        theta1_approxCI_up <- sapply(theta1_approxCI_up[1], Rmpfr::asNumeric)})
        #theta1_approxCI_up <- Rmpfr::asNumeric(theta1_approxCI_up)})
        
      }else{
        
        T1 <- data[data<tau[1]]
        n1 <- length(T1)
        #T2 <- data[data>=tau[1]] #Avner: Not used later so not necessary (also would be problematic considering non-censored data)
        
        p1 <- 1 - exp(-tau[1]/mle1)
        p2 <- exp(-tau[1]/mle1) * (1 - p1*exp(-(tau[2]-tau[1])/theta21)  - (1-p1)*exp(-(tau[2]-tau[1])/theta22))
              exp(-tau[1]/mle1) * (1 - p1*exp(-(tau[2]-tau[1])/theta21)  - (1-p1)*exp(-(tau[2]-tau[1])/theta22))
        p3 <- 1 - p1 - p2
        Cn <- 1/(1 - (1-p1)^n - (1-p2)^n + p3^n)
        
        # add p1, p3
        bias_part <- function(n){
          y <- 0
          for (i in 1:(n-1)) {
            for (k in 0:i) {
              C_ik <- (-1)^k * choose(n, i) * choose(i, k) * ((1 - p1)^(n-i) - p3^(n-i)) * (1-p1)^k
              tau_ik <- tau[1] / i * (n - i + k)
              re_ik <- C_ik *tau_ik
              y <- y + re_ik
            }
          }
          return(y)
        }
        
        V1 <- mle1^2/n1
        
        theta1_approxCI_low <- mle1 - Cn * bias_part(n) - qnorm(1-alpha/2) * sqrt(V1)
        theta1_approxCI_low <- ifelse(theta1_approxCI_low < 0, 0L, theta1_approxCI_low)
        theta1_approxCI_up <- mle1 - Cn * bias_part(n) + qnorm(1-alpha/2) * sqrt(V1)
      }
      
      
      ################# Theta21,22 p Type I
      
      
      # sample_dat_2nd is the observations under the second stress level
      # n1, n2 are the number of failures under two stress levels respectively
      # p, theta21, theta22 are the estimates from the EM algorithm
      # Compute the elements in the hessian matrix
      
      sample_dat_2nd <- data[data>=tau[1]]
      sample_dat_2nd <- sample_dat_2nd[sample_dat_2nd < tau[2]]
      
      n2 <- length(sample_dat_2nd)
      n_f <- length(data[data<tau[2]])
      
      loglik_i_cont <- function(x){
        theta21 <- x[1]
        theta22 <- x[2]
        p <- x[3]
        #d is the data
        sum(log(p/theta21*exp(-(sample_dat_2nd-tau[1])/theta21)+(1-p)/theta22*exp(-(sample_dat_2nd-tau[1])/theta22))) +
          (n-n_f)*log(p*exp(-(tau[2]-tau[1])/theta21)+(1-p)*exp(-(tau[2]-tau[1])/theta22))
      }

      x <- c(theta21,theta22,p)
      hess <- numDeriv::hessian(func = loglik_i_cont, x = x)
      
      
    }else{######## Type II ###########
      
      
      if(n>25){
        #####################################
        ###The case with a relatively large sample size (e.g., n = 50)
        ###High precision computation is needed
        #T1_simu <- mpfr(rexp(n, rate = 1/theta1), 256) # 256 can be increased
        suppressWarnings({
        T1_simu <- Rmpfr::mpfr(1, 256)*data
        n1 <- sum(T1_simu <= tau)
        T1 <- sort(T1_simu[T1_simu <= tau])
        mle1 <- Rmpfr::mpfr(sum(T1, (n-n1)*tau) / n1, 256)
        ###The modified version of bias
        bias_Mpfr_modify <- function(x, n, r){
          y <- 0
          p_j_sum <- 0
          for (i in 1:(r-1)) {
            med <- Rmpfr::chooseMpfr(n, i)*(1 - exp(-tau/x))^i*exp(-tau/x)^(n-i)
            p_j_sum <- p_j_sum + med
          }
          for (j in 1:(r-1)) {
            for (k in 0:j) {
              c_jk <- (-1)^k / p_j_sum * Rmpfr::chooseMpfr(n, j) * Rmpfr::chooseMpfr(j, k) * exp(-tau*(n - j + k)/x)
              tau_jk <- tau / j * (n - j + k)
              re_jk <- c_jk *tau_jk
              y <- y + re_jk
            }
          }
          return(y)
        }
        
        V1 <- mle1^2/n1
        
        bias <- bias_Mpfr_modify(mle1, n = Rmpfr::mpfr(n, 256), r = r)  #n = mpfr(n, 512) when n = 200
        theta1_approxCI_low <- mle1 - bias - qnorm(1-alpha/2)*sqrt(V1)
        theta1_approxCI_low <- ifelse(theta1_approxCI_low < 0, 0L, theta1_approxCI_low)
        #Avner: Reverted to older version due to error message
        theta1_approxCI_low <- sapply(theta1_approxCI_low[1], Rmpfr::asNumeric)
        #theta1_approxCI_low <- Rmpfr::asNumeric(theta1_approxCI_low)
        theta1_approxCI_up <- mle1 - bias + qnorm(1-alpha/2)*sqrt(V1)
        #Avner: Same as above
        theta1_approxCI_up <- sapply(theta1_approxCI_up[1], Rmpfr::asNumeric)})
        #theta1_approxCI_up <- Rmpfr::asNumeric(theta1_approxCI_up)})
      }else{
        
        T1 <- data[data<tau[1]]
        n1 <- length(T1)
        
        ####theta1
        bias <- function(x, r){
          y <- NULL
          val1 <- 1 - exp(-tau/x)
          p_j_sum <- pbinom(0, n, val1, lower.tail = FALSE) - pbinom(r-1, n, val1, lower.tail = FALSE)
          for (j in 1:(r-1)) {
            for (k in 0:j) {
              c_jk <- (-1)^k / p_j_sum * choose(n, j) * choose(j, k) * exp(-tau*(n - j + k)/x)
              tau_jk <- tau / j * (n - j + k)
              re_jk <- c_jk *tau_jk
              y <- c(y, re_jk)
            }
          }
          return(sum(y))
        }
        
        V1 <- mle1^2/n1
        
        theta1_approxCI_low <- mle1-bias(mle1, r) - qnorm(1-alpha/2)*sqrt(V1)
        theta1_approxCI_low <- ifelse(theta1_approxCI_low < 0, 0L, theta1_approxCI_low)
        theta1_approxCI_up <- mle1-bias(mle1, r) + qnorm(1-alpha/2)*sqrt(V1)
        
      }
      
      
      ################# Theta21,22 p Type II
      
      # sample_dat_2nd is the observations under the second stress level
      sample_dat_2nd <- data[data>=tau[1]]
      # tr is the rth ordered failure, we expect data to be ordered
      tr <- sample_dat_2nd[r-n1] #Avner: Doesn't seem to be used after its defined.
      sample_dat_2nd <- sample_dat_2nd[1:(r-n1)]
      
      loglik_ii_cont <- function(x){
        theta21 <- x[1]
        theta22 <- x[2]
        p <- x[3]
        sum(log(p/theta21*exp(-(sample_dat_2nd-tau)/theta21)+(1-p)/theta22*exp(-(sample_dat_2nd-tau)/theta22))) +
          (n-r)*log(p*exp(-(max(sample_dat_2nd)-tau)/theta21)+(1-p)*exp(-(max(sample_dat_2nd)-tau)/theta22))
      }
      
      x <- c(theta21,theta22,p)
      hess <- numDeriv::hessian(func = loglik_ii_cont, x = x)
      
    }
    
    
  }
  
  ###Get the inverse, the diagonal elements are the variances of the parameters 
  CV_matrice <- solve(-hess)
  Vtheta21 <- diag(CV_matrice)[1] 
  Vtheta22 <- diag(CV_matrice)[2] 
  Vp <- diag(CV_matrice)[3]
  
  ### alpha is the significance level
  
  p_approxCI_low <- p - qnorm(1-alpha/2)*sqrt(Vp)
  p_approxCI_low <- ifelse(p_approxCI_low < 0, 0L, p_approxCI_low) #compare it to 0 and 1!
  p_approxCI_up <- p + qnorm(1-alpha/2)*sqrt(Vp)
  p_approxCI_up <- ifelse(p_approxCI_up > 1, 1L, p_approxCI_up)
  
  if (n>25 && monitoring=="continuous") {
    p_approxCI_low <- sapply(p_approxCI_low[1], Rmpfr::asNumeric)
    p_approxCI_up <- sapply(p_approxCI_up[1], Rmpfr::asNumeric)
  }
  
  theta21_approxCI_low <- theta21 - qnorm(1-alpha/2)*sqrt(Vtheta21)
  theta21_approxCI_low <- ifelse(theta21_approxCI_low < 0, 0L, theta21_approxCI_low)
  theta21_approxCI_up <- theta21 + qnorm(1-alpha/2)*sqrt(Vtheta21)
  
  theta22_approxCI_low <- theta22 - qnorm(1-alpha/2)*sqrt(Vtheta22)
  theta22_approxCI_low <- ifelse(theta22_approxCI_low < 0, 0L, theta22_approxCI_low)
  theta22_approxCI_up <- theta22 + qnorm(1-alpha/2)*sqrt(Vtheta22)
  
  #Confidence Intervals
  theta1_CI <- c(theta1_approxCI_low,theta1_approxCI_up)
  theta21_CI <- c(theta21_approxCI_low,theta21_approxCI_up)
  theta22_CI <- c(theta22_approxCI_low,theta22_approxCI_up)
  p_CI <- c(p_approxCI_low,p_approxCI_up)
  
  output <- list(theta1=theta1_CI, theta21=theta21_CI, theta22=theta22_CI,p=p_CI, B=NA, j=NA, alpha=alpha, type="Asymptotic")
  class(output) <- "CIhSSALT"
  return(output)
  
  
}