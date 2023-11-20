

CIsay_hSSALT <- function(data, n, censoring, tau , r, monitoring, delta, alpha, theta1, p, theta21, theta22){
  
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
    
    if (n > sum(data)) { #data is censored with a tau2 that is smaller than the input tau2
      n2j <- data[-(1:q1)] # censored
      if(q2 > length(n2j)){n2j <- c(n2j, rep(0, q2-length(n2j)))}
    }else{n2j <- data[-c((1:q1),((q1+q2+1):length(data)))]}
    
    n2 <- sum(n2j)
    
    
    A <- A_func(p = p, theta21 = theta21, theta22 = theta22, tau2j0, tau2j, tau1 = tau[1])
    B <- B_func(p = p, theta21 = theta21, theta22 = theta22, tau1 = tau[1], tau2 = tau[2])
    
    nf <- n1+n2
    O22 <- -Secondderivative_O22(p , theta21, theta22 , tau2j0 = tau2j0, tau2j = tau2j, n2j = n2j, nf = nf, tau1 = tau[1], tau2 = tau[2], A, B)
    O23 <- -Secondderivative_O23(p , theta21, theta22 , tau2j0 = tau2j0, tau2j = tau2j, n2j = n2j, nf = nf, tau1 = tau[1], tau2 = tau[2], A, B)
    O24 <- -Secondderivative_O24(p , theta21, theta22 , tau2j0 = tau2j0, tau2j = tau2j, n2j = n2j, nf = nf, tau1 = tau[1], tau2 = tau[2], A, B)
    O33 <- -Secondderivative_O33(p , theta21, theta22 , tau2j0 = tau2j0, tau2j = tau2j, n2j = n2j, nf = nf, tau1 = tau[1], tau2 = tau[2], A, B)
    O34 <- -Secondderivative_O34(p , theta21, theta22 , tau2j0 = tau2j0, tau2j = tau2j, n2j = n2j, nf = nf, tau1 = tau[1], tau2 = tau[2], A, B)
    O44 <- -Secondderivative_O44(p , theta21, theta22 , tau2j0 = tau2j0, tau2j = tau2j, n2j = n2j, nf = nf, tau1 = tau[1], tau2 = tau[2], A, B)
    
    
    O11 <- -Secondderivative_O11(theta1 = theta1, tau1j, tau1j0, n1j, n1, tau1 = tau[1])
    V_theta1 <- 1/O11
    V_theta1
    
    # alpha is the significant level 
    theta1_approxCI_low <- theta1 - qnorm(1-alpha/2)*sqrt(V_theta1)
    # theta1_approxCI_low #compare it to 0
    theta1_approxCI_up <- theta1 + qnorm(1-alpha/2)*sqrt(V_theta1)
    # theta1_approxCI_up
    
    
    ####Compute the inverse of the Fisher matrix and get the covariance matrix
    CV_matrice <- solve(matrix(c(O22, O23, O24, O23, O33, O34, O24, O34, O44), nrow = 3, byrow = TRUE))
    V_theta21 <- diag(CV_matrice)[1] 
    V_theta22 <- diag(CV_matrice)[2] 
    V_p <- diag(CV_matrice)[3] 
    
    theta21_approxCI_low <- theta21 - qnorm(1-alpha/2)*sqrt(V_theta21)
    # theta21_approxCI_low #compare it to 0
    theta21_approxCI_up <- theta21 + qnorm(1-alpha/2)*sqrt(V_theta21)

    
    theta22_approxCI_low <- theta22 - qnorm(1-alpha/2)*sqrt(V_theta22)
    # theta22_approxCI_low #compare it to 0
    theta22_approxCI_up <- theta22 + qnorm(1-alpha/2)*sqrt(V_theta22)
    # theta22_approxCI_up
    
    p_approxCI_low <- p - qnorm(1-alpha/2)*sqrt(V_p)
    # p_approxCI_low #compare it to 0 and 1!
    p_approxCI_up <- p + qnorm(1-alpha/2)*sqrt(V_p)
    # p_approxCI_up
    
    return(list(CI_theta_1 = c(theta1_approxCI_low, theta1_approxCI_up),
                   CI_theta_21 = c(theta21_approxCI_low, theta21_approxCI_up), 
                   CI_theta_22 = c(theta22_approxCI_low, theta22_approxCI_up),
                   CI_theta_p = c(p_approxCI_low, p_approxCI_up)
                   ))
    
    
    
    
    
  }else{############ continuous ###############
    
  
    mle1 <- theta1
    
    
    
    
    if(censoring==1){
      
      # Theta1, Type I
      
      if(n>50){
        
        ###The case with a relatively large sample size (e.g., n = 50)
        ###High precision computation is needed
        library(Rmpfr)
        set.seed(1234)  ### The data should also be newly simulated!!! Otherwise, different precision level causes problem
        T1_simu <- mpfr(rexp(n, rate = 1/theta1), 256) # 64 can be increased
        n1 <- sum(T1_simu <= tau[1])
        T1 <- sort(T1_simu[T1_simu <= tau[1]])
        mle1 <- mpfr(sum(T1, (n-n1)*tau[1]) / n1, 256)
        
        # add mle1, tau
        bias_part_Mpfr <- function(n){
          y <- 0
          p1 <- 1 - exp(-tau[1]/mle1)
          p2 <- exp(-tau[1]/mle1) * (1 - p*exp(-(tau[2]-tau[1])/theta21) - (1-p)*exp(-(tau[2]-tau[1])/theta22))
          p3 <- 1 - p1 - p2
          for (i in 1:(n-1)) {
            for (k in 0:i) {
              n <- mpfr(n, 256)
              C_ik <- (-1)^k * chooseMpfr(n, i) * chooseMpfr(i, k) * ((1 - p1)^(n-i) - p3^(n-i)) * (1-p1)^k
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
        
        # alpha <- c(0.1, 0.05, 0.01)
        theta1_approxCI_low <- mle1 - Cn * bias - qnorm(1-alpha/2)*sqrt(V1)
        theta1_approxCI_low <- ifelse(theta1_approxCI_low < 0, 0L, theta1_approxCI_low)
        theta1_approxCI_up <- mle1 - Cn * bias + qnorm(1-alpha/2)*sqrt(V1)
        
        
      }else{
        
        T1 <- data[data<tau[1]]
        n1 <- length(T1)
        T2 <- data[data>=tau[1]]
        
        
        
        p1 <- 1 - exp(-tau[1]/mle1)
        p2 <- exp(-tau[1]/mle1) * (1 - p1*exp(-(tau[2]-tau[1])/theta21)  - (1-p1)*exp(-(tau[2]-tau[1])/theta22))
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
        # V1
        
        # alpha <- c(0.1, 0.05, 0.01)
        theta1_approxCI_low <- mle1 - Cn * bias_part(n) - qnorm(1-alpha/2) * sqrt(V1)
        theta1_approxCI_low <- ifelse(theta1_approxCI_low < 0, 0L, theta1_approxCI_low)
        # theta1_approxCI_low 
        theta1_approxCI_up <- mle1 - Cn * bias_part(n) + qnorm(1-alpha/2) * sqrt(V1)
        # theta1_approxCI_up
      }
      
      
      
      
      
      ################# Theta21,22 p Type I
      
      
      # sample_dat_2nd is the observations under the second stress level
      # n1, n2 are the number of failures under two stress levels respectively
      # p, theta21, theta22 are the estimates from the EM algorithm
      # Compute the elements in the hessian matrix
      
      
      sample_dat_2nd <- data[data>=tau[1]]
      sample_dat_2nd <- sample_dat_2nd[sample_dat_2nd < tau[2]]
      
      n2 <- length(sample_dat_2nd)
      
      Otheta21 <- O_theta21_I(sample_dat_2nd, p, theta21, theta22, n1, n2, tau1 = tau[1], tau2 = tau[2])
      Otheta21theta22 <- O_theta21theta22_I(sample_dat_2nd, p, theta21, theta22, n1, n2, tau1 = tau[1], tau2 = tau[2])
      Otheta21p <- O_theta21p_I(sample_dat_2nd, p, theta21, theta22, n1, n2, tau1 = tau[1], tau2 = tau[2])
      Otheta22 <- O_theta22_I(sample_dat_2nd, p, theta21, theta22, n1, n2, tau1 = tau[1], tau2 = tau[2])
      Otheta22p <- O_theta22p_I(sample_dat_2nd, p, theta21, theta22, n1, n2, tau1 = tau[1], tau2 = tau[2])
      Op <- O_p_I(sample_dat_2nd, p, theta21, theta22, n1, n2, tau1 = tau[1], tau2 = tau[2])
      
      ###Get the inverse, the diagonal elements are the variances of the parameters 
      CV_matrice <- solve(matrix(c(Otheta21, Otheta21theta22, Otheta21p, Otheta21theta22, 
                                   Otheta22, Otheta22p, Otheta21p, Otheta22p, Op), nrow = 3, byrow = TRUE))
      Vtheta21 <- diag(CV_matrice)[1] 
      Vtheta22 <- diag(CV_matrice)[2] 
      Vp <- diag(CV_matrice)[3] 
      
      ### alpha is the significant level
      p_approxCI_low <- p1 - qnorm(1-alpha/2)*sqrt(Vp)
      p_approxCI_low <- ifelse(p_approxCI_low < 0, 0L, p_approxCI_low) #compare it to 0 and 1!
      # p_approxCI_low
      p_approxCI_up <- p1 + qnorm(1-alpha/2)*sqrt(Vp)
      # p_approxCI_up
      
      
      theta21_approxCI_low <- theta21 - qnorm(1-alpha/2)*sqrt(Vtheta21)
      theta21_approxCI_low <- ifelse(theta21_approxCI_low < 0, 0L, theta21_approxCI_low)
      # theta21_approxCI_low
      theta21_approxCI_up <- theta21 + qnorm(1-alpha/2)*sqrt(Vtheta21)
      # theta21_approxCI_up
      
      theta22_approxCI_low <- theta22 - qnorm(1-alpha/2)*sqrt(Vtheta22)
      theta22_approxCI_low <- ifelse(theta22_approxCI_low < 0, 0L, theta22_approxCI_low)
      # theta22_approxCI_low
      theta22_approxCI_up <- theta22 + qnorm(1-alpha/2)*sqrt(Vtheta22)
      # theta22_approxCI_up
     
      
      
      
      
      
      
       
    }else{######## Type II ###########
      
      
      if(n>50){
        #####################################
        ###The case with a relatively large sample size (e.g., n = 50)
        ###High precision computation is needed
        library(Rmpfr)
        set.seed(12345)   ### The data should also be newly simulated!!! Otherwise, different precision level causes problem
        T1_simu <- mpfr(rexp(n, rate = 1/theta1), 256) # 256 can be increased
        n1 <- sum(T1_simu <= tau)
        T1 <- sort(T1_simu[T1_simu <= tau])
        mle1 <- mpfr(sum(T1, (n-n1)*tau) / n1, 256)
        ###The modified version of bias
        bias_Mpfr_modify <- function(x, n, r){
          y <- 0
          p_j_sum <- 0
          for (i in 1:(r-1)) {
            med <- chooseMpfr(n, i)*(1 - exp(-tau/x))^i*exp(-tau/x)^(n-i)
            p_j_sum <- p_j_sum + med
          }
          for (j in 1:(r-1)) {
            for (k in 0:j) {
              c_jk <- (-1)^k / p_j_sum * chooseMpfr(n, j) * chooseMpfr(j, k) * exp(-tau*(n - j + k)/x)
              tau_jk <- tau / j * (n - j + k)
              re_jk <- c_jk *tau_jk
              y <- y + re_jk
            }
          }
          return(y)
        }
        
        V1 <- mle1^2/n1
        V1
        
        bias <- bias_Mpfr_modify(mle1, n = mpfr(n, 256), r = r)  #n = mpfr(n, 512) when n = 200
        # alpha <- c(0.1, 0.05, 0.01)
        theta1_approxCI_low <- mle1 - bias - qnorm(1-alpha/2)*sqrt(V1)
        theta1_approxCI_low <- ifelse(theta1_approxCI_low < 0, 0L, theta1_approxCI_low)
        # theta1_approxCI_low  
        theta1_approxCI_up <- mle1 - bias + qnorm(1-alpha/2)*sqrt(V1)
        # theta1_approxCI_up
      }else{
        
        T1 <- data[data<tau[1]]
        n1 <- length(T1)
        # T2 <- data[data>=tau[1]]
      
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
      # V1
      
      # alpha <- c(0.1, 0.05, 0.01)
      theta1_approxCI_low <- mle1-bias(mle1, r) - qnorm(1-alpha/2)*sqrt(V1)
      theta1_approxCI_low <- ifelse(theta1_approxCI_low < 0, 0L, theta1_approxCI_low)
      # theta1_approxCI_low 
      theta1_approxCI_up <- mle1-bias(mle1, r) + qnorm(1-alpha/2)*sqrt(V1)
      # theta1_approxCI_up
      

      }
      
      
      ################# Theta21,22 p Type II
      
      # sample_dat_2nd is the observations under the second stress level
      sample_dat_2nd <- data[data>=tau[1]]
      # tr is the rth ordered failure, we expect data to be ordered
      tr <- sample_dat_2nd[r-n1]
      sample_dat_2nd <- sample_dat_2nd[1:(r-n1)]

      
      # p, theta21, theta22 are the estimates from the EM algorithm, inputted
      
      
      
      # Compute the elements in the hessien matrix
      
      
      Otheta21 <- O_theta21_II(ti = sample_dat_2nd, p = p, theta21 = theta21, 
                            theta22 = theta22, tr = tr, r = r, tau)
      Otheta21theta22 <- O_theta21theta22_II(ti = sample_dat_2nd, p = p, theta21 = theta21, 
                                          theta22 = theta22, tr = tr, r = r, tau)
      Otheta21p <- O_theta21p_II(ti = sample_dat_2nd, p = p, theta21 = theta21, 
                              theta22 = theta22, tr = tr, r = r, tau)
      Otheta22 <- O_theta22_II(ti = sample_dat_2nd, p = p, theta21 = theta21, 
                            theta22 = theta22, tr = tr, r = r, tau)
      Otheta22p <-  O_theta22p_II(ti = sample_dat_2nd, p = p, theta21 = theta21, 
                               theta22 = theta22, tr = tr, r = r, tau)
      Op <- O_p_II(ti = sample_dat_2nd, p = p, theta21 = theta21, 
                theta22 = theta22, tr = tr, r = r, tau)
      
      ###Get the inverse, the diagonal elements are the variances of the parameters 
      CV_matrice <- solve(matrix(c(Otheta21, Otheta21theta22, Otheta21p, Otheta21theta22, 
                                   Otheta22, Otheta22p, Otheta21p, Otheta22p, Op), nrow = 3, byrow = TRUE))
      Vtheta21 <- diag(CV_matrice)[1] 
      Vtheta22 <- diag(CV_matrice)[2] 
      Vp <- diag(CV_matrice)[3] 
      
      ### alpha is the significant level
      p_approxCI_low <- p - qnorm(1-alpha/2)*sqrt(Vp)
      p_approxCI_low <- ifelse(p_approxCI_low < 0, 0L, p_approxCI_low) #compare it to 0 and 1!
      p_approxCI_low
      p_approxCI_up <- p + qnorm(1-alpha/2)*sqrt(Vp)
      p_approxCI_up
      
      
      theta21_approxCI_low <- theta21 - qnorm(1-alpha/2)*sqrt(Vtheta21)
      theta21_approxCI_low <- ifelse(theta21_approxCI_low < 0, 0L, theta21_approxCI_low)
      theta21_approxCI_low
      theta21_approxCI_up <- theta21 + qnorm(1-alpha/2)*sqrt(Vtheta21)
      theta21_approxCI_up
      
      theta22_approxCI_low <- theta22 - qnorm(1-alpha/2)*sqrt(Vtheta22)
      theta22_approxCI_low <- ifelse(theta22_approxCI_low < 0, 0L, theta22_approxCI_low)
      theta22_approxCI_low
      theta22_approxCI_up <- theta22 + qnorm(1-alpha/2)*sqrt(Vtheta22)
      theta22_approxCI_up
      
    }
    
    
    return(list(CI_theta_1 = c(theta1_approxCI_low, theta1_approxCI_up),
                   CI_theta_21 = c(theta21_approxCI_low, theta21_approxCI_up), 
                   CI_theta_22 = c(theta22_approxCI_low, theta22_approxCI_up),
                   CI_theta_p = c(p_approxCI_low, p_approxCI_up)
    ))
    
    
  }
  
  
  
  
  
  
}