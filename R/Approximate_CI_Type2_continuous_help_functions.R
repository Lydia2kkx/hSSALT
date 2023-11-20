
####Approximate CIs
####Continuous Case
####Type-II censoring 
####Parameter under the second stress level
#######################################################################
################ 3 parameters in the mixture
#######################################################################
# O_theta21 is the negative second derivative of theta21
# A1, A2, B1, B2 are four parts to compute O_theta21
A1_II <- function(p, theta21, theta22, ti, tau){
  nominator <- p/(theta21^3)*exp(-(ti-tau)/theta21)*(((ti-tau)/theta21)^2 - 4*(ti-tau)/theta21 + 2)
  denominator <- p/theta21*exp(-(ti-tau)/theta21) + (1-p)/theta22*exp(-(ti-tau)/theta22)
  y <- sum(nominator/denominator)
  return(y)
}

A2_II <- function(p, theta21, theta22, ti, tau){
  nominator <- p/(theta21^2)*exp(-(ti-tau)/theta21)*((ti-tau)/theta21 -1)
  denominator <- p/theta21*exp(-(ti-tau)/theta21) + (1-p)/theta22*exp(-(ti-tau)/theta22)
  y <- sum((nominator/denominator)^2)
  return(y)
}

B1_II <- function(p, theta21, theta22, tr, r, tau){
  nominator <- (n-r)*(tr-tau)*p/(theta21^3)*exp(-(tr-tau)/theta21)*((tr-tau)/theta21 -2)
  denominator <- p*exp(-(tr-tau)/theta21) + (1-p)*exp(-(tr-tau)/theta22)
  y <- nominator/denominator
  return(y)
}

B2_II <- function(p, theta21, theta22, tr, r, tau){
  nominator <- (tr-tau)*p/(theta21^2)*exp(-(tr-tau)/theta21)
  denominator <- p*exp(-(tr-tau)/theta21) + (1-p)*exp(-(tr-tau)/theta22)
  y <- (n-r)*(nominator/denominator)^2
  return(y)
}

O_theta21_II <- function(ti, p, theta21, theta22, tr, r, tau){
  O22 <- -(A1_II(p, theta21, theta22, ti, tau) - A2_II(p, theta21, theta22, ti, tau) + 
             B1_II(p, theta21, theta22, tr, r, tau) - B2_II(p, theta21, theta22, tr, r, tau))
  return(O22)
}


####################################################################################
# O_theta21theta22 is the negative second derivative of theta21 and theta22
# A23_1, B23_1 are two parts to compute O_theta21theta22
#######################O23
A23_1_II <- function(p, theta21, theta22, ti, tau){
  nominator <- p/(theta21^2)*(1-p)/(theta22^2)*exp(-(ti-tau)/theta21)*exp(-(ti-tau)/theta22)*((ti-tau)/theta21-1)*((ti-tau)/theta22-1)
  denominator <- p/theta21*exp(-(ti-tau)/theta21) + (1-p)/theta22*exp(-(ti-tau)/theta22)
  y <- sum(nominator/(denominator^2))
  return(y)
}

B23_1_II <- function(p, theta21, theta22, tr, r, tau){
  nominator <- (n-r)*(tr-tau)^2*p/(theta21^2)*(1-p)/(theta22^2)*exp(-(tr-tau)/theta21)*exp(-(tr-tau)/theta22)
  denominator <- p*exp(-(tr-tau)/theta21) + (1-p)*exp(-(tr-tau)/theta22)
  y <- nominator/(denominator^2)
  return(y)
}

O_theta21theta22_II <- function(ti, p, theta21, theta22, tr, r, tau){
  O23 <- A23_1_II(p, theta21, theta22, ti, tau) + B23_1_II(p, theta21, theta22, tr, r, tau) 
  return(O23)
}

####################################################################################
# O_theta21p is the negative second derivative of theta21 and p
# A24_1, A24_2, B24_1, B24_2 are four parts to compute O_theta21p
#######################O24
A24_1_II <- function(p, theta21, theta22, ti, tau){
  nominator <- 1/(theta21^2)*exp(-(ti-tau)/theta21)*((ti-tau)/theta21-1)
  denominator <- p/theta21*exp(-(ti-tau)/theta21) + (1-p)/theta22*exp(-(ti-tau)/theta22)
  y <- sum(nominator/denominator)
  return(y)
}

A24_2_II <- function(p, theta21, theta22, ti, tau){
  nominator <- p/(theta21^2)*exp(-(ti-tau)/theta21)*((ti-tau)/theta21-1)*(exp(-(ti-tau)/theta21)/theta21 - exp(-(ti-tau)/theta22)/theta22)
  denominator <- p/theta21*exp(-(ti-tau)/theta21) + (1-p)/theta22*exp(-(ti-tau)/theta22)
  y <- sum(nominator/(denominator^2))
  return(y)
}

B24_1_II <- function(p, theta21, theta22, tr, r, tau){
  nominator <- (n-r)*(tr-tau)/(theta21^2)*exp(-(tr-tau)/theta21)
  denominator <- p*exp(-(tr-tau)/theta21) + (1-p)*exp(-(tr-tau)/theta22)
  y <- nominator/denominator
  return(y)
}

B24_2_II <- function(p, theta21, theta22, tr, r, tau){
  nominator <- (n-r)*(tr-tau)*p/(theta21^2)*exp(-(tr-tau)/theta21)*(exp(-(tr-tau)/theta21)-exp(-(tr-tau)/theta22))
  denominator <- p*exp(-(tr-tau)/theta21) + (1-p)*exp(-(tr-tau)/theta22)
  y <- nominator/(denominator^2)
  return(y)
}

O_theta21p_II <- function(ti, p, theta21, theta22, tr, r, tau){
  O24 <- -A24_1_II(p, theta21, theta22, ti, tau) + A24_2_II(p, theta21, theta22, ti, tau) - 
    B24_1_II(p, theta21, theta22, tr, r, tau) + B24_2_II(p, theta21, theta22, tr, r, tau) 
  return(O24)
}

####################################################################################
# O_theta22 is the negative second derivative of theta22
# A1_2, A2_2, B1_2, B2_2 are four parts to compute O_theta22
#######################O33
####For theta22
A1_2_II <- function(p, theta21, theta22, ti, tau){
  nominator <- (1-p)/(theta22^3)*exp(-(ti-tau)/theta22)*(((ti-tau)/theta22)^2 - 4*(ti-tau)/theta22 + 2)
  denominator <- p/theta21*exp(-(ti-tau)/theta21) + (1-p)/theta22*exp(-(ti-tau)/theta22)
  y <- sum(nominator/denominator)
  return(y)
}

A2_2_II <- function(p, theta21, theta22, ti, tau){
  nominator <- (1-p)/(theta22^2)*exp(-(ti-tau)/theta22)*((ti-tau)/theta22 -1)
  denominator <- p/theta21*exp(-(ti-tau)/theta21) + (1-p)/theta22*exp(-(ti-tau)/theta22)
  y <- sum((nominator/denominator)^2)
  return(y)
}

B1_2_II <- function(p, theta21, theta22, tr, r, tau){
  nominator <- (n-r)*(tr-tau)*(1-p)/(theta22^3)*exp(-(tr-tau)/theta22)*((tr-tau)/theta22 -2)
  denominator <- p*exp(-(tr-tau)/theta21) + (1-p)*exp(-(tr-tau)/theta22)
  y <- nominator/denominator
  return(y)
}

B2_2_II <- function(p, theta21, theta22, tr, r, tau){
  nominator <- (tr-tau)*(1-p)/(theta22^2)*exp(-(tr-tau)/theta22)
  denominator <- p*exp(-(tr-tau)/theta21) + (1-p)*exp(-(tr-tau)/theta22)
  y <- (n-r)*(nominator/denominator)^2
  return(y)
}

O_theta22_II <- function(ti, p, theta21, theta22, tr, r, tau){
  O33 <- -(A1_2_II(p, theta21, theta22, ti, tau) - A2_2_II(p, theta21, theta22, ti, tau) + 
             B1_2_II(p, theta21, theta22, tr, r, tau) - B2_2_II(p, theta21, theta22, tr, r, tau))
  return(O33)
}

####################################################################################
# O_theta22p is the negative second derivative of theta22 and p
# A34_1, A34_2, B34_1, B34_2 are four parts to compute O_theta22p
#######################O34
A34_1_II <- function(p, theta21, theta22, ti, tau){
  nominator <- 1/(theta22^2)*exp(-(ti-tau)/theta22)*((ti-tau)/theta22-1)
  denominator <- p/theta21*exp(-(ti-tau)/theta21) + (1-p)/theta22*exp(-(ti-tau)/theta22)
  y <- sum(nominator/denominator)
  return(y)
}

A34_2_II <- function(p, theta21, theta22, ti, tau){
  nominator <- (1-p)/(theta22^2)*exp(-(ti-tau)/theta22)*((ti-tau)/theta22-1)*(exp(-(ti-tau)/theta21)/theta21 - exp(-(ti-tau)/theta22)/theta22)
  denominator <- p/theta21*exp(-(ti-tau)/theta21) + (1-p)/theta22*exp(-(ti-tau)/theta22)
  y <- sum(nominator/(denominator^2))
  return(y)
}

B34_1_II <- function(p, theta21, theta22, tr, r, tau){
  nominator <- (n-r)*(tr-tau)/(theta22^2)*exp(-(tr-tau)/theta22)
  denominator <- p*exp(-(tr-tau)/theta21) + (1-p)*exp(-(tr-tau)/theta22)
  y <- nominator/denominator
  return(y)
}

B34_2_II <- function(p, theta21, theta22, tr, r, tau){
  nominator <- (n-r)*(tr-tau)*(1-p)/(theta22^2)*exp(-(tr-tau)/theta22)*(exp(-(tr-tau)/theta21)-exp(-(tr-tau)/theta22))
  denominator <- p*exp(-(tr-tau)/theta21) + (1-p)*exp(-(tr-tau)/theta22)
  y <- nominator/(denominator^2)
  return(y)
}

O_theta22p_II <- function(ti, p, theta21, theta22, tr, r, tau){
  O34 <- A34_1_II(p, theta21, theta22, ti, tau) + A34_2_II(p, theta21, theta22, ti, tau) + 
    B34_1_II(p, theta21, theta22, tr, r, tau) + B34_2_II(p, theta21, theta22, tr, r, tau) 
  return(O34)
}

####################################################################################
# O_p is the negative second derivative of p
# A_p, B_p, B34_1, B34_2 are two parts to compute O_p
######For p
A_p_II <- function(p, theta21, theta22, ti, tau){
  nominator <- 1/theta21*exp(-(ti-tau)/theta21) - 1/theta22*exp(-(ti-tau)/theta22)
  denominator <- p/theta21*exp(-(ti-tau)/theta21) + (1-p)/theta22*exp(-(ti-tau)/theta22)
  y <- sum((nominator/denominator)^2)
  return(y)
}

B_p_II <- function(p, theta21, theta22, tr, r, tau){
  nominator <- exp(-(tr-tau)/theta21) - exp(-(tr-tau)/theta22)
  denominator <- p*exp(-(tr-tau)/theta21) + (1-p)*exp(-(tr-tau)/theta22)
  y <- (n-r)*(nominator/denominator)^2
  return(y)
}

O_p_II <- function(ti, p, theta21, theta22, tr, r, tau){
  O44 <- A_p_II(p, theta21, theta22, ti, tau) + B_p_II(p, theta21, theta22, tr, r, tau)
  return(O44)
}