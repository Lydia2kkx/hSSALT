####Approximate CIs
####Interval Case
######Formulas for Fisher Information Matrix for interval monitoring (have sent you)
######2 stress levels are computed together, no need to do bias adjustment
######Oij means the position of the element

Secondderivative_O11 <- function(theta1, tau1j, tau1j0, n1j, n1, tau1){
  -2/(theta1^3)*sum(n1j*(tau1j0*exp(-tau1j0/theta1)-tau1j*exp(-tau1j/theta1))/(exp(-tau1j0/theta1)-exp(-tau1j/theta1))) + 
    1/(theta1^2)*sum(n1j*((tau1j0/theta1)^2*exp(-tau1j0/theta1)-(tau1j/theta1)^2*exp(-tau1j/theta1))/(exp(-tau1j0/theta1)-exp(-tau1j/theta1)))-
    1/(theta1^2)*sum(n1j*((tau1j0/theta1*exp(-tau1j0/theta1)-tau1j/theta1*exp(-tau1j/theta1))/(exp(-tau1j0/theta1)-exp(-tau1j/theta1)))^2)-
    2*tau1*(n-n1)/(theta1^3)
}

A_func <- function(p, theta21, theta22, tau2j0, tau2j, tau1){
  p*(exp(-(tau2j0-tau1)/theta21)-exp(-(tau2j-tau1)/theta21)) + 
    (1-p)*(exp(-(tau2j0-tau1)/theta22)-exp(-(tau2j-tau1)/theta22))
}

B_func <- function(p, theta21, theta22, tau1, tau2){
  p*exp(-(tau2-tau1)/theta21) + (1-p)*exp(-(tau2-tau1)/theta22)
}

Secondderivative_O22 <- function(p, theta21, theta22, tau2j0, tau2j, n2j, nf, tau1, tau2, A, B){
  -2*p/(theta21^3)*sum(n2j*((tau2j0-tau1)*exp(-(tau2j0-tau1)/theta21)-(tau2j-tau1)*exp(-(tau2j-tau1)/theta21))/A) + 
    p/(theta21^2)*sum(n2j*((tau2j0-tau1)^2*exp(-(tau2j0-tau1)/theta21)-(tau2j-tau1)^2*exp(-(tau2j-tau1)/theta21))/(A*theta21^2))-
    p/(theta21^2)*sum(n2j*p*(((tau2j0-tau1)/theta21*exp(-(tau2j0-tau1)/theta21)-(tau2j-tau1)/theta21*exp(-(tau2j-tau1)/theta21))/A)^2)-
    2*p/(theta21^3)*(n-nf)*(tau2-tau1)*exp(-(tau2-tau1)/theta21)/B + 
    p/(theta21^2)*(n-nf)*((tau2-tau1)/theta21)^2*exp(-(tau2-tau1)/theta21)/B - 
    (p/theta21)^2*(n-nf)*((tau2-tau1)/theta21*exp(-(tau2-tau1)/theta21)/B)^2
}


Secondderivative_O33 <- function(p, theta21, theta22, tau2j0, tau2j, n2j, nf, tau1, tau2, A, B){
  -2*(1-p)/(theta22^3)*sum(n2j*((tau2j0-tau1)*exp(-(tau2j0-tau1)/theta22)-(tau2j-tau1)*exp(-(tau2j-tau1)/theta22))/A) + 
    (1-p)/(theta22^2)*sum(n2j*((tau2j0-tau1)^2*exp(-(tau2j0-tau1)/theta22)-(tau2j-tau1)^2*exp(-(tau2j-tau1)/theta22))/(A*theta22^2))-
    (1-p)/(theta22^2)*sum(n2j*(1-p)*(((tau2j0-tau1)/theta22*exp(-(tau2j0-tau1)/theta22)-(tau2j-tau1)/theta22*exp(-(tau2j-tau1)/theta22))/A)^2)-
    2*(1-p)/(theta22^3)*(n-nf)*(tau2-tau1)*exp(-(tau2-tau1)/theta22)/B + 
    (1-p)/(theta22^2)*(n-nf)*((tau2-tau1)/theta22)^2*exp(-(tau2-tau1)/theta22)/B - 
    ((1-p)/theta22)^2*(n-nf)*((tau2-tau1)/theta22*exp(-(tau2-tau1)/theta22)/B)^2
}


Secondderivative_O44 <- function(p, theta21, theta22, tau2j0, tau2j, n2j, nf, tau1, tau2, A, B){
  -sum(n2j*(((exp(-(tau2j0-tau1)/theta21)- exp(-(tau2j-tau1)/theta21))-(exp(-(tau2j0-tau1)/theta22)- exp(-(tau2j-tau1)/theta22)))/A)^2) -
    (n-nf)*((exp(-(tau2-tau1)/theta21)-exp(-(tau2-tau1)/theta22))/B)^2
}


Secondderivative_O23 <- function(p, theta21, theta22, tau2j0, tau2j, n2j, nf, tau1, tau2, A, B){
  -sum(n2j*(p/(theta21^2)*((tau2j0-tau1)*exp(-(tau2j0-tau1)/theta21)-(tau2j-tau1)*exp(-(tau2j-tau1)/theta21))/A)*
         ((1-p)/(theta22^2)*((tau2j0-tau1)*exp(-(tau2j0-tau1)/theta22)-(tau2j-tau1)*exp(-(tau2j-tau1)/theta22))/A))-
    (n-nf)*(tau2-tau1)^2*p/(theta21^2)*exp(-(tau2-tau1)/theta21)*(1-p)/(theta22^2)*exp(-(tau2-tau1)/theta22)/(B^2)
}

Secondderivative_O24 <- function(p, theta21, theta22, tau2j0, tau2j, n2j, nf, tau1, tau2, A, B){
  sum(n2j/(theta21^2)*((tau2j0-tau1)*exp(-(tau2j0-tau1)/theta21)-(tau2j-tau1)*exp(-(tau2j-tau1)/theta21))/A)-
    sum(n2j*p/(theta21^2)*((tau2j0-tau1)*exp(-(tau2j0-tau1)/theta21)-(tau2j-tau1)*exp(-(tau2j-tau1)/theta21))/A*
          ((exp(-(tau2j0-tau1)/theta21)- exp(-(tau2j-tau1)/theta21))-(exp(-(tau2j0-tau1)/theta22)- exp(-(tau2j-tau1)/theta22)))/A)+
    (n-nf)*(tau2-tau1)/(theta21^2)*exp(-(tau2-tau1)/theta21)/B-
    (n-nf)*(tau2-tau1)*p/(theta21^2)*exp(-(tau2-tau1)/theta21)*(exp(-(tau2-tau1)/theta21)-exp(-(tau2-tau1)/theta22))/(B^2)
}

Secondderivative_O34 <- function(p, theta21, theta22, tau2j0, tau2j, n2j, nf, tau1, tau2, A, B){
  -sum(n2j/(theta22^2)*((tau2j0-tau1)*exp(-(tau2j0-tau1)/theta22)-(tau2j-tau1)*exp(-(tau2j-tau1)/theta22))/A)-
    sum(n2j*(1-p)/(theta22^2)*((tau2j0-tau1)*exp(-(tau2j0-tau1)/theta22)-(tau2j-tau1)*exp(-(tau2j-tau1)/theta22))/A*
          ((exp(-(tau2j0-tau1)/theta21)- exp(-(tau2j-tau1)/theta21))-(exp(-(tau2j0-tau1)/theta22)- exp(-(tau2j-tau1)/theta22)))/A)-
    (n-nf)*(tau2-tau1)/(theta22^2)*exp(-(tau2-tau1)/theta22)/B-
    (n-nf)*(tau2-tau1)*(1-p)/(theta22^2)*exp(-(tau2-tau1)/theta22)*(exp(-(tau2-tau1)/theta21)-exp(-(tau2-tau1)/theta22))/(B^2)
}