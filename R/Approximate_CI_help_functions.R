loglik_i_interval <- function(x){
  theta21 <- x[1]
  theta22 <- x[2]
  p <- x[3]
  
  value1 <- p*(exp(-(tau2j0-tau[1])/theta21)-exp(-(tau2j-tau[1])/theta21))
  value2 <- (1-p)*(exp(-(tau2j0-tau[1])/theta22)-exp(-(tau2j-tau[1])/theta22))
  
  sum(n2j*log(value1+value2)) +
    (n-nf)*log(p*exp(-(tau[2]-tau[1])/theta21)+(1-p)*exp(-(tau[2]-tau[1])/theta22))
  
}

loglik_i_cont <- function(x){
  theta21 <- x[1]
  theta22 <- x[2]
  p <- x[3]
  #d is the data
  sum(log(p/theta21*exp(-(sample_dat_2nd-tau[1])/theta21)+(1-p)/theta22*exp(-(sample_dat_2nd-tau[1])/theta22))) +
    (n-n_f)*log(p*exp(-(tau[2]-tau[1])/theta21)+(1-p)*exp(-(tau[2]-tau[1])/theta22))
}

loglik_ii_cont <- function(x){
  theta21 <- x[1]
  theta22 <- x[2]
  p <- x[3]
  sum(log(p/theta21*exp(-(sample_dat_2nd-tau)/theta21)+(1-p)/theta22*exp(-(sample_dat_2nd-tau)/theta22))) +
    (n-r)*log(p*exp(-(max(sample_dat_2nd)-tau)/theta21)+(1-p)*exp(-(max(sample_dat_2nd)-tau)/theta22))
}

Secondderivative_O11 <- function(theta1, tau1j, tau1j0, n1j, n1, tau1){
  -2/(theta1^3)*sum(n1j*(tau1j0*exp(-tau1j0/theta1)-tau1j*exp(-tau1j/theta1))/(exp(-tau1j0/theta1)-exp(-tau1j/theta1))) + 
    1/(theta1^2)*sum(n1j*((tau1j0/theta1)^2*exp(-tau1j0/theta1)-(tau1j/theta1)^2*exp(-tau1j/theta1))/(exp(-tau1j0/theta1)-exp(-tau1j/theta1)))-
    1/(theta1^2)*sum(n1j*((tau1j0/theta1*exp(-tau1j0/theta1)-tau1j/theta1*exp(-tau1j/theta1))/(exp(-tau1j0/theta1)-exp(-tau1j/theta1)))^2)-
    2*tau1*(n-n1)/(theta1^3)
}