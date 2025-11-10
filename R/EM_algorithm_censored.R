
###Find the result with the largest log likelihood
find_max <- function(x){
  index <- which.max(x$loglik)
  y <- x[index, ]
  row.names(y) <- NULL
  return(y)
}

###Sum up the finite elements
sum_finite <- function(x) {
  sum(x[is.finite(x)])
}

###EM algorithm for continuous monitoring 
EM_algorithm_censored <- function(ind, data, d, N, parameter_starts, tol){

  p1 <- parameter_starts[ind, 1] 
  p2 <-  1 - p1
  rate1 <- 1/parameter_starts[ind, 2]
  rate2 <- 1/parameter_starts[ind, 3]

  loglik<- rep(NA, N)
  loglik[1]<-0
  loglik[2]<-sum_finite(p1*(log(p1) + d*log(dexp(data, rate1)) + (1-d)*log(1-pexp(data, rate1))))+
    sum_finite(p2*(log(p2) + d*log(dexp(data, rate2))+ (1-d)*log(1-pexp(data, rate2))))
  
  k<-2
  
  while((abs(loglik[k]-loglik[k-1]) >= tol) && (k <= N)) {
    # E step
    sum.of.comps1 <- p1*dexp(data, rate1)+p2*dexp(data, rate2)
    sum.of.comps2 <- p1*(1-pexp(data, rate1))+p2*(1-pexp(data, rate2))
    hij1 <- (p1*dexp(data, rate1)/sum.of.comps1)^d*(p1*(1-pexp(data, rate1))/sum.of.comps2)^(1-d)
    hij2 <- (p2*dexp(data, rate2)/sum.of.comps1)^d*(p2*(1-pexp(data, rate2))/sum.of.comps2)^(1-d)
    
    # M step
    p1<-sum_finite(hij1)/length(data)
    p2<-sum_finite(hij2)/length(data)
    
    rate1 <- sum_finite(hij1*d)/sum_finite(hij1*data)
    rate2 <- sum_finite(hij2*d)/sum_finite(hij2*data)
    
    loglik[k+1]<- sum_finite(hij1*(log(p1)+d*log(dexp(data, rate1)) + (1-d)*log(1-pexp(data, rate1))))+
      sum_finite(hij2*(log(p2)+d*log(dexp(data, rate2))+ (1-d)*log(1-pexp(data, rate2))))
    k<-k+1
  }
  
  ###theta21 should be smaller than theta22 in the output
  p_helper <- p1 
  if(rate1 <= rate2){
    theta21 <- 1/rate2
    theta22 <- 1/rate1
    p1 <- p2
    p2 <- p_helper
  }else{
    theta21 <- 1/rate1
    theta22 <- 1/rate2
  }
  
  if(k <= N){
    message <- "convergent"
    posterior <- cbind(hij1, hij2)
  } else {
    message <- "not convergent"
    posterior <- NA
    warning("The EM Algorithm failed to converge!")
  }
  
  return(list(results = data.frame(p1 = p1, p2 = p2, theta21 = theta21, theta22 = theta22,
                                     loglik = loglik[k], iteration = k-1,
                                     message = message), posterior = posterior))
}




