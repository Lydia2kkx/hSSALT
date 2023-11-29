# misc


### Find the result with the largest log likelihood
find_max <- function(x){
  index <- which.max(x$loglik)
  y <- x[index, ]
  row.names(y) <- NULL
  return(y)
}


mysum <- function(x) {
  sum(x[is.finite(x)])
}



### Compute mle1 in interval case
func_theta1_int <- function(x, n ,n1j, n1, tau1, tau1j, tau1j0){
  ### Vector calculation
  sum(n1j*(exp(-tau1j0/x)*tau1j0 - exp(-tau1j/x)*tau1j) / (exp(-tau1j0/x) - exp(-tau1j/x))) + tau1*(n-n1)
}






EM_algorithm_censored <- function(ind, data, d, N, parameter_starts, tol){
  pi1 <- parameter_starts[ind, 1]
  pi2 <-  1 - pi1
  rate1 <- 1/parameter_starts[ind, 2]
  rate2 <- 1/parameter_starts[ind, 3]

  # print("outside")
  # print(c(pi1,pi2,rate1,rate2))



  loglik<- rep(NA, N)
  loglik[1]<-0
  loglik[2]<-mysum(pi1*(log(pi1) + d*log(dexp(data, rate1)) + (1-d)*log(1-pexp(data, rate1))))+
    mysum(pi2*(log(pi2) + d*log(dexp(data, rate2))+ (1-d)*log(1-pexp(data, rate2))))

  k<-2

  # loop
  while((abs(loglik[k]-loglik[k-1]) >= tol) && (k <= N)) {
    # E step
    sum.of.comps1 <- pi1*dexp(data, rate1)+pi2*dexp(data, rate2)
    sum.of.comps2 <- pi1*(1-pexp(data, rate1))+pi2*(1-pexp(data, rate2))
    tau1 <- (pi1*dexp(data, rate1)/sum.of.comps1)^d*(pi1*(1-pexp(data, rate1))/sum.of.comps2)^(1-d)
    tau2 <- (pi2*dexp(data, rate2)/sum.of.comps1)^d*(pi2*(1-pexp(data, rate2))/sum.of.comps2)^(1-d)

    # M step
    pi1<-mysum(tau1)/length(data)
    pi2<-mysum(tau2)/length(data)

    rate1 <- mysum(tau1*d)/mysum(tau1*data)
    rate2 <- mysum(tau2*d)/mysum(tau2*data)

    loglik[k+1]<- mysum(tau1*(log(pi1)+d*log(dexp(data, rate1)) + (1-d)*log(1-pexp(data, rate1))))+
      mysum(tau2*(log(pi2)+d*log(dexp(data, rate2))+ (1-d)*log(1-pexp(data, rate2))))
    k<-k+1
  }
  if(k <= N){
    p_helper <- pi1
    if(rate1 <= rate2){
      theta21 <- 1/rate2
      theta22 <- 1/rate1
      pi1 <- pi2
      pi2 <- p_helper
      posterior <- cbind(tau2, tau1)
    }else{
      theta21 <- 1/rate1
      theta22 <- 1/rate2
      posterior <- cbind(tau1, tau2)
    }
    return(list(results = data.frame(pi1 = pi1, pi2 = pi2, theta21 = theta21, theta22 = theta22,
                      loglik = loglik[k], iteration = k-1,
                      message = "convergent"), posterior = posterior))
  }else{
    return(list(results = data.frame(pi1 = NA, pi2 = NA, theta21 = NA, theta22 = NA,
                      loglik = NA, iteration = k-1,
                      message = "not convergent"), posterior = NA))
  }
}



EM_algorithm_interval <- function(ind, data , N, delta , d, parameter_starts, q2, tol){
  omega1 <- parameter_starts[ind, 1]
  omega2 <-  1 - omega1
  p1 <- 1 - exp(-delta/parameter_starts[ind, 2])
  p2 <- 1 - exp(-delta/parameter_starts[ind, 3])
  mysum <- function(x) {
    sum(x[is.finite(x)])
  }
  loglik<- rep(NA, N)
  loglik[1]<-0
  loglik[2]<-mysum(omega1*(log(omega1)+d*log(dgeom(data, p1))+(1-d)*log(1-pgeom(data, p1))))+
    mysum(omega2*(log(omega2)+log(dgeom(data, p2))+(1-d)*log(1-pgeom(data, p2))))

  k<-2

  # loop
  while((abs(loglik[k]-loglik[k-1]) >= tol) & (k <= N)) {
    # E step
    sum.of.comps1 <- omega1*dgeom(data, p1)+omega2*dgeom(data, p2)
    sum.of.comps2 <- omega1*(1-pgeom(data, p1))+omega2*(1-pgeom(data, p2))
    tau1 <- (omega1*dgeom(data, p1)/sum.of.comps1)^d*(omega1*(1-pgeom(data, p1))/sum.of.comps2)^(1-d)
    tau2 <- (omega2*dgeom(data, p2)/sum.of.comps1)^d*(omega2*(1-pgeom(data, p2))/sum.of.comps2)^(1-d)

    # M step
    omega1 <- mysum(tau1)/mysum(tau1 + tau2)
    omega2 <- mysum(tau2)/mysum(tau1 + tau2)

    p1 <- mysum(tau1*d)/mysum(tau1*(d*(data+1) + (1-d)*q2))
    p2 <- mysum(tau2*d)/mysum(tau2*(d*(data+1) + (1-d)*q2))

    loglik[k+1]<-mysum(tau1*(log(omega1)+d*log(dgeom(data, p1)) + (1-d)*log(1-pgeom(data, p1))))+
      mysum(tau2*(log(omega2)+d*log(dgeom(data, p2))+ (1-d)*log(1-pgeom(data, p2))))

    k<-k+1
  }
  if(k <= N){
    p_helper <- omega1
    if(p1 <= p2){
      theta21 <- -delta/log(1 - p2)
      theta22 <- -delta/log(1 - p1)
      omega1 <- omega2
      omega2 <- p_helper
      posterior <- cbind(tau2, tau1)
    }else{
      theta21 <- -delta/log(1 - p1)
      theta22 <- -delta/log(1 - p2)
      posterior <- cbind(tau1, tau2)
    }
    return(list(results = data.frame(ind = ind, prob1 = omega1, prob2 = omega2, theta21 = theta21, theta22 = theta22,
                      loglik = loglik[k], iteration = k-1, message = "convergent"),
                posterior = posterior))
  }else{
    return(list(results = data.frame(ind = ind, prob1 = NA, prob2 = NA, theta21 = NA, theta22 = NA,
                      loglik = NA, iteration = k-1, message = "not convergent"),
                posterior = NA))
  }
}
