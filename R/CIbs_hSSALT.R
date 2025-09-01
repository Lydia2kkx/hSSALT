
CIbs_hSSALT <-function(data, n, censoring, tau , r, monitoring, delta, alpha, B, theta1, theta21,
                       theta22, p, maxit, tol, language, parallel, ncores,grid){

    
  bootstrap_distri_list <- bootstrap_distribution(data, n, monitoring = monitoring, theta1 = theta1,
                                             theta21 = theta21, theta22 = theta22, p = p, 
                                             censoring, tau, r, B, delta = delta, maxit, tol, 
                                             language, parallel, ncores, grid)
  
  
  bootstrap_distri <- bootstrap_distri_list[[1]]
  num_of_iterations <- bootstrap_distri_list[[2]]
  
  
  theta1_low <- quantile(bootstrap_distri$theta1, alpha/2)[[1]]
  theta1_up <- quantile(bootstrap_distri$theta1, 1-alpha/2)[[1]]
  
  theta21_low <- quantile(bootstrap_distri$theta21, alpha/2)[[1]]
  theta21_up <- quantile(bootstrap_distri$theta21, 1-alpha/2)[[1]]
    
  theta22_low <- quantile(bootstrap_distri$theta22, alpha/2)[[1]]
  theta22_up <- quantile(bootstrap_distri$theta22, 1-alpha/2)[[1]]
  
  p_low <- quantile(bootstrap_distri$p1, alpha/2)[[1]]
  p_up <- quantile(bootstrap_distri$p1, 1-alpha/2)[[1]]
  
  #Confidence Intervals
  theta1_CI <- c(theta1_low,theta1_up)
  theta21_CI <- c(theta21_low,theta21_up)
  theta22_CI <- c(theta22_low,theta22_up)
  p_CI <- c(p_low,p_up)
  
  output <- list(theta1=theta1_CI, theta21=theta21_CI, theta22=theta22_CI,p=p_CI, B=B, j=num_of_iterations, alpha=alpha, type="Percentile")
  class(output) <- "CIhSSALT"
  return(output)

}
