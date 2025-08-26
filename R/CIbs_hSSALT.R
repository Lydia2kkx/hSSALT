
CIbs_hSSALT <-function(data, n, censoring, tau , r, monitoring, delta, alpha, B, theta1, theta21,
                       theta22, p, maxit, tol, language, parallel, ncores,grid){

    
  bootstrap_distri <- bootstrap_distribution(data, n, monitoring = monitoring, theta1 = theta1,
                                             theta21 = theta21, theta22 = theta22, p = p, 
                                             censoring, tau, r, B, delta = delta, maxit, tol, 
                                             language, parallel, ncores, grid)
    
  
  lower_bound <- c(quantile(bootstrap_distri$theta1, alpha/2)[[1]],quantile(bootstrap_distri$theta21, alpha/2)[[1]],quantile(bootstrap_distri$theta22, alpha/2)[[1]],quantile(bootstrap_distri$p1, alpha/2)[[1]])
  upper_bound <- c(quantile(bootstrap_distri$theta1, 1-alpha/2)[[1]],quantile(bootstrap_distri$theta21, 1-alpha/2)[[1]],quantile(bootstrap_distri$theta22, 1-alpha/2)[[1]],quantile(bootstrap_distri$p1, 1-alpha/2)[[1]])
  
  conf_ints <- cbind(lower_bound,upper_bound)
  row.names(conf_ints) <- c("theta1","theta21","theta22","p")
  
  return(conf_ints)

}
