### Compute mle1 in interval case
func_theta1_int <- function(x, n ,n1j, n1, hijk1, tau1j, tau1j0){
  ### Vector calculation
  sum(n1j*(exp(-tau1j0/x)*tau1j0 - exp(-tau1j/x)*tau1j) / (exp(-tau1j0/x) - exp(-tau1j/x))) + hijk1*(n-n1)
}

EM_algorithm_interval <- function(ind, data , N, delta , d, parameter_starts, q2, tol){
  
  omega1 <- parameter_starts[ind, 1]
  omega2 <-  1 - omega1
  p1 <- 1 - exp(-delta/parameter_starts[ind, 2])
  p2 <- 1 - exp(-delta/parameter_starts[ind, 3])
  sum_finite <- function(x) {
    sum(x[is.finite(x)])
  }
  loglik<- rep(NA, N)
  loglik[1]<-0
  loglik[2]<-sum_finite(omega1*(log(omega1)+d*log(dgeom(data, p1))+(1-d)*log(1-pgeom(data, p1))))+
    sum_finite(omega2*(log(omega2)+log(dgeom(data, p2))+(1-d)*log(1-pgeom(data, p2))))
  
  k<-2
  
  # loop
  while((abs(loglik[k]-loglik[k-1]) >= tol) & (k <= N)) {
    # E step
    sum.of.comps1 <- omega1*dgeom(data, p1)+omega2*dgeom(data, p2)
    sum.of.comps2 <- omega1*(1-pgeom(data, p1))+omega2*(1-pgeom(data, p2))
    hijk1 <- (omega1*dgeom(data, p1)/sum.of.comps1)^d*(omega1*(1-pgeom(data, p1))/sum.of.comps2)^(1-d)
    hijk2 <- (omega2*dgeom(data, p2)/sum.of.comps1)^d*(omega2*(1-pgeom(data, p2))/sum.of.comps2)^(1-d)
    
    # M step
    omega1 <- sum_finite(hijk1)/sum_finite(hijk1 + hijk2)
    omega2 <- sum_finite(hijk2)/sum_finite(hijk1 + hijk2)
    
    p1 <- sum_finite(hijk1*d)/sum_finite(hijk1*(d*(data+1) + (1-d)*q2))
    p2 <- sum_finite(hijk2*d)/sum_finite(hijk2*(d*(data+1) + (1-d)*q2))
    
    loglik[k+1]<-sum_finite(hijk1*(log(omega1)+d*log(dgeom(data, p1)) + (1-d)*log(1-pgeom(data, p1))))+
      sum_finite(hijk2*(log(omega2)+d*log(dgeom(data, p2))+ (1-d)*log(1-pgeom(data, p2))))
    
    k<-k+1
  }
  
  p_helper <- omega1
  if(p1 <= p2){
    theta21 <- -delta/log(1 - p2)
    theta22 <- -delta/log(1 - p1)
    omega1 <- omega2
    omega2 <- p_helper
    posterior <- cbind(hijk2, hijk1)
  }else{
    theta21 <- -delta/log(1 - p1)
    theta22 <- -delta/log(1 - p2)
    posterior <- cbind(hijk1, hijk2)
  }
  
  if(k <= N){
    message <- "convergent"
  } else {
    message <- "not convergent"
    posterior <- NA
    warning("The EM Algorithm failed to converge!")
  }
  
  return(list(results = data.frame(p1 = omega1, p2 = omega2, theta21 = theta21, theta22 = theta22,
                                   loglik = loglik[k], iteration = k-1,
                                   message = message), posterior = posterior))
}