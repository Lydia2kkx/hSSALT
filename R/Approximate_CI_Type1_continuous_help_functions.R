####Approximate CIs
####Continuous Case
####Type-I censoring 
####Parameter under the second stress level
#######################################################################
################ 3 parameters in the mixture
#######################################################################
# O_theta21 is the negative second derivative of theta21
# A1, A2, B1, B2 are four parts to compute O_theta21
A1_I <- function(p, theta21, theta22, ti, tau1){
  nominator <- p/(theta21^3)*exp(-(ti-tau1)/theta21)*(((ti-tau1)/theta21)^2 - 4*(ti-tau1)/theta21 + 2)
  denominator <- p/theta21*exp(-(ti-tau1)/theta21) + (1-p)/theta22*exp(-(ti-tau1)/theta22)
  y <- sum(nominator/denominator)
  return(y)
}

A2_I <- function(p, theta21, theta22, ti, tau1){
  nominator <- p/(theta21^2)*exp(-(ti-tau1)/theta21)*((ti-tau1)/theta21 -1)
  denominator <- p/theta21*exp(-(ti-tau1)/theta21) + (1-p)/theta22*exp(-(ti-tau1)/theta22)
  y <- sum((nominator/denominator)^2)
  return(y)
}

B1_I <- function(p, theta21, theta22, n1, n2, tau1, tau2){
  nominator <- (n-n1-n2)*(tau2-tau1)*p/(theta21^3)*exp(-(tau2-tau1)/theta21)*((tau2-tau1)/theta21 -2)
  denominator <- p*exp(-(tau2-tau1)/theta21) + (1-p)*exp(-(tau2-tau1)/theta22)
  y <- nominator/denominator
  return(y)
}

B2_I <- function(p, theta21, theta22, n1, n2, tau1, tau2){
  nominator <- (tau2-tau1)*p/(theta21^2)*exp(-(tau2-tau1)/theta21)
  denominator <- p*exp(-(tau2-tau1)/theta21) + (1-p)*exp(-(tau2-tau1)/theta22)
  y <- (n-n1-n2)*(nominator/denominator)^2
  return(y)
}

O_theta21_I <- function(ti, p, theta21, theta22, n1, n2, tau1, tau2){
  O22 <- -(A1_I(p, theta21, theta22, ti, tau1) - A2_I(p, theta21, theta22, ti, tau1) + 
             B1_I(p, theta21, theta22, n1, n2, tau1, tau2) - B2_I(p, theta21, theta22, n1, n2,tau1, tau2))
  return(O22)
}

####################################################################################
# O_theta21theta22 is the negative second derivative of theta21 and theta22
# A23_1, B23_1 are two parts to compute O_theta21theta22
#######################O23
A23_1_I <- function(p, theta21, theta22, ti, tau1){
  nominator <- p/(theta21^2)*(1-p)/(theta22^2)*exp(-(ti-tau1)/theta21)*exp(-(ti-tau1)/theta22)*((ti-tau1)/theta21-1)*((ti-tau1)/theta22-1)
  denominator <- p/theta21*exp(-(ti-tau1)/theta21) + (1-p)/theta22*exp(-(ti-tau1)/theta22)
  y <- sum(nominator/(denominator^2))
  return(y)
}

B23_1_I <- function(p, theta21, theta22, n1, n2, tau1, tau2){
  nominator <- (n-n1-n2)*(tau2-tau1)^2*p/(theta21^2)*(1-p)/(theta22^2)*exp(-(tau2-tau1)/theta21)*exp(-(tau2-tau1)/theta22)
  denominator <- p*exp(-(tau2-tau1)/theta21) + (1-p)*exp(-(tau2-tau1)/theta22)
  y <- nominator/(denominator^2)
  return(y)
}

O_theta21theta22_I <- function(ti, p, theta21, theta22, n1, n2, tau1, tau2){
  O23 <- A23_1_I(p, theta21, theta22, ti, tau1) + B23_1_I(p, theta21, theta22, n1, n2, tau1, tau2) 
  return(O23)
}

####################################################################################
# O_theta21p is the negative second derivative of theta21 and p
# A24_1, A24_2, B24_1, B24_2 are four parts to compute O_theta21p
#######################O24
A24_1_I <- function(p, theta21, theta22, ti, tau1){
  nominator <- 1/(theta21^2)*exp(-(ti-tau1)/theta21)*((ti-tau1)/theta21-1)
  denominator <- p/theta21*exp(-(ti-tau1)/theta21) + (1-p)/theta22*exp(-(ti-tau1)/theta22)
  y <- sum(nominator/denominator)
  return(y)
}

A24_2_I <- function(p, theta21, theta22, ti, tau1){
  nominator <- p/(theta21^2)*exp(-(ti-tau1)/theta21)*((ti-tau1)/theta21-1)*(exp(-(ti-tau1)/theta21)/theta21 - exp(-(ti-tau1)/theta22)/theta22)
  denominator <- p/theta21*exp(-(ti-tau1)/theta21) + (1-p)/theta22*exp(-(ti-tau1)/theta22)
  y <- sum(nominator/(denominator^2))
  return(y)
}

B24_1_I <- function(p, theta21, theta22,n1, n2, tau1, tau2){
  nominator <- (n-n1-n2)*(tau2-tau1)/(theta21^2)*exp(-(tau2-tau1)/theta21)
  denominator <- p*exp(-(tau2-tau1)/theta21) + (1-p)*exp(-(tau2-tau1)/theta22)
  y <- nominator/denominator
  return(y)
}

B24_2_I <- function(p, theta21, theta22,n1, n2, tau1, tau2){
  nominator <- (n-n1-n2)*(tau2-tau1)*p/(theta21^2)*exp(-(tau2-tau1)/theta21)*(exp(-(tau2-tau1)/theta21)-exp(-(tau2-tau1)/theta22))
  denominator <- p*exp(-(tau2-tau1)/theta21) + (1-p)*exp(-(tau2-tau1)/theta22)
  y <- nominator/(denominator^2)
  return(y)
}

O_theta21p_I <- function(ti, p, theta21, theta22,n1, n2, tau1, tau2){
  O24 <- -A24_1_I(p, theta21, theta22, ti, tau1) + A24_2_I(p, theta21, theta22, ti, tau1) - 
    B24_1_I(p, theta21, theta22,n1, n2, tau1, tau2) + B24_2_I(p, theta21, theta22,n1, n2, tau1, tau2) 
  return(O24)
}

####################################################################################
# O_theta22 is the negative second derivative of theta22
# A1_2, A2_2, B1_2, B2_2 are four parts to compute O_theta22
#######################O33
####For theta22
A1_2_I <- function(p, theta21, theta22, ti, tau1){
  nominator <- (1-p)/(theta22^3)*exp(-(ti-tau1)/theta22)*(((ti-tau1)/theta22)^2 - 4*(ti-tau1)/theta22 + 2)
  denominator <- p/theta21*exp(-(ti-tau1)/theta21) + (1-p)/theta22*exp(-(ti-tau1)/theta22)
  y <- sum(nominator/denominator)
  return(y)
}

A2_2_I <- function(p, theta21, theta22, ti, tau1){
  nominator <- (1-p)/(theta22^2)*exp(-(ti-tau1)/theta22)*((ti-tau1)/theta22 -1)
  denominator <- p/theta21*exp(-(ti-tau1)/theta21) + (1-p)/theta22*exp(-(ti-tau1)/theta22)
  y <- sum((nominator/denominator)^2)
  return(y)
}

B1_2_I <- function(p, theta21, theta22,n1, n2, tau1, tau2){
  nominator <- (n-n1-n2)*(tau2-tau1)*(1-p)/(theta22^3)*exp(-(tau2-tau1)/theta22)*((tau2-tau1)/theta22 -2)
  denominator <- p*exp(-(tau2-tau1)/theta21) + (1-p)*exp(-(tau2-tau1)/theta22)
  y <- nominator/denominator
  return(y)
}

B2_2_I <- function(p, theta21, theta22,n1, n2, tau1, tau2){
  nominator <- (tau2-tau1)*(1-p)/(theta22^2)*exp(-(tau2-tau1)/theta22)
  denominator <- p*exp(-(tau2-tau1)/theta21) + (1-p)*exp(-(tau2-tau1)/theta22)
  y <- (n-n1-n2)*(nominator/denominator)^2
  return(y)
}

O_theta22_I <- function(ti, p, theta21, theta22,n1, n2, tau1, tau2){
  O33 <- -(A1_2_I(p, theta21, theta22, ti, tau1) - A2_2_I(p, theta21, theta22, ti, tau1) + 
             B1_2_I(p, theta21, theta22,n1, n2, tau1, tau2) - B2_2_I(p, theta21, theta22,n1, n2, tau1, tau2))
  return(O33)
}

####################################################################################
# O_theta22p is the negative second derivative of theta22 and p
# A34_1, A34_2, B34_1, B34_2 are four parts to compute O_theta22p
#######################O34
A34_1_I <- function(p, theta21, theta22, ti, tau1){
  nominator <- 1/(theta22^2)*exp(-(ti-tau1)/theta22)*((ti-tau1)/theta22-1)
  denominator <- p/theta21*exp(-(ti-tau1)/theta21) + (1-p)/theta22*exp(-(ti-tau1)/theta22)
  y <- sum(nominator/denominator)
  return(y)
}

A34_2_I <- function(p, theta21, theta22, ti, tau1){
  nominator <- (1-p)/(theta22^2)*exp(-(ti-tau1)/theta22)*((ti-tau1)/theta22-1)*(exp(-(ti-tau1)/theta21)/theta21 - exp(-(ti-tau1)/theta22)/theta22)
  denominator <- p/theta21*exp(-(ti-tau1)/theta21) + (1-p)/theta22*exp(-(ti-tau1)/theta22)
  y <- sum(nominator/(denominator^2))
  return(y)
}

B34_1_I <- function(p, theta21, theta22,n1, n2, tau1, tau2){
  nominator <- (n-n1-n2)*(tau2-tau1)/(theta22^2)*exp(-(tau2-tau1)/theta22)
  denominator <- p*exp(-(tau2-tau1)/theta21) + (1-p)*exp(-(tau2-tau1)/theta22)
  y <- nominator/denominator
  return(y)
}

B34_2_I <- function(p, theta21, theta22,n1, n2, tau1, tau2){
  nominator <- (n-n1-n2)*(tau2-tau1)*(1-p)/(theta22^2)*exp(-(tau2-tau1)/theta22)*(exp(-(tau2-tau1)/theta21)-exp(-(tau2-tau1)/theta22))
  denominator <- p*exp(-(tau2-tau1)/theta21) + (1-p)*exp(-(tau2-tau1)/theta22)
  y <- nominator/(denominator^2)
  return(y)
}

O_theta22p_I <- function(ti, p, theta21, theta22,n1, n2, tau1, tau2){
  O34 <- A34_1_I(p, theta21, theta22, ti, tau1) + A34_2_I(p, theta21, theta22, ti, tau1) + 
    B34_1_I(p, theta21, theta22,n1, n2, tau1, tau2) + B34_2_I(p, theta21, theta22,n1, n2, tau1, tau2) 
  return(O34)
}

####################################################################################
# O_p is the negative second derivative of p
# A_p, B_p, B34_1, B34_2 are two parts to compute O_p
######For p
A_p_I <- function(p, theta21, theta22, ti, tau1){
  nominator <- 1/theta21*exp(-(ti-tau1)/theta21) - 1/theta22*exp(-(ti-tau1)/theta22)
  denominator <- p/theta21*exp(-(ti-tau1)/theta21) + (1-p)/theta22*exp(-(ti-tau1)/theta22)
  y <- sum((nominator/denominator)^2)
  return(y)
}

B_p_I <- function(p, theta21, theta22,n1, n2, tau1, tau2){
  nominator <- exp(-(tau2-tau1)/theta21) - exp(-(tau2-tau1)/theta22)
  denominator <- p*exp(-(tau2-tau1)/theta21) + (1-p)*exp(-(tau2-tau1)/theta22)
  y <- (n-n1-n2)*(nominator/denominator)^2
  return(y)
}

O_p_I <- function(ti, p, theta21, theta22,n1, n2, tau1, tau2){
  O44 <- A_p_I(p, theta21, theta22, ti, tau1) + B_p_I(p, theta21, theta22,n1, n2, tau1, tau2)
  return(O44)
}