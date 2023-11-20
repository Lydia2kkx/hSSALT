#include "em_algorithm_int.h"



Rcpp::DataFrame EM_algorithm_interval_arma(const arma::vec& data, const double& delta, const double& ind, const arma::vec& d, const Rcpp::DataFrame& parameter_starts, const double& q2, const int& N, const double& tol){
  arma::vec firstcol = parameter_starts[0];
  arma::vec secondcol = parameter_starts[1];
  arma::vec thirdcol = parameter_starts[2];
  
  double omega1 = firstcol[ind-1];
  
  double omega2 =  1 - omega1;
  double p1 = 1 - exp(-delta/secondcol[(ind-1)]);
  double p2 = 1 - exp(-delta/thirdcol[(ind-1)]);
  
  
  
  arma::vec loglik(N);
  loglik[0] = 0;
  
  loglik[1]=mysum_cpp(omega1 * (log(omega1)+d%log(dgeom_boost(data, p1))+(1-d)%log(1-pgeom_boost(data, p1))))+ 
    mysum_cpp(omega2 * (log(omega2)+log(dgeom_boost(data, p2))+(1-d)%log(1-pgeom_boost(data, p2))));          
  
  int k =1;
  
  // loop
  while ((abs(loglik[k]-loglik[k-1]) >= tol) && (k <= (N-1))) {
   
    
    // E step
    arma::vec sum_of_comps1 = omega1 * dgeom_boost(data, p1)+omega2 * dgeom_boost(data, p2); 
    arma::vec sum_of_comps2 = omega1 * (1-pgeom_boost(data, p1))+omega2 * (1-pgeom_boost(data, p2));  
    arma::vec tau1 = arma::pow((omega1 * dgeom_boost(data, p1)/sum_of_comps1),d)%arma::pow((omega1 * (1-pgeom_boost(data, p1))/sum_of_comps2),(1-d)); 
    arma::vec tau2 = arma::pow((omega2 * dgeom_boost(data, p2)/sum_of_comps1),d)%arma::pow((omega2 * (1-pgeom_boost(data, p2))/sum_of_comps2),(1-d)); 
    
    
    
    // M step
    omega1 = mysum_cpp(tau1)/mysum_cpp(tau1 + tau2); 
    omega2 = mysum_cpp(tau2)/mysum_cpp(tau1 + tau2); 
    
    p1 = mysum_cpp(tau1%d)/mysum_cpp(tau1%(d%(data+1) + (1-d) * q2));
    p2 = mysum_cpp(tau2%d)/mysum_cpp(tau2%(d%(data+1) + (1-d) * q2));
    
   
    
    loglik[k+1]=mysum_cpp(tau1%(log(omega1)+d%arma::log(dgeom_boost(data, p1)) + (1-d)%arma::log(1-pgeom_boost(data, p1))))+ 
      mysum_cpp(tau2%(log(omega2)+d%arma::log(dgeom_boost(data, p2))+ (1-d)%arma::log(1-pgeom_boost(data, p2))));           
    
    
    
    k++;
    
    
  }
  if(k <= (N-1)){
    double p_helper = omega1;
    if(p1 <= p2){
      double theta21 = -delta/std::log(1 - p2);
      double theta22 = -delta/std::log(1 - p1);
      double omega1 = omega2;
      double omega2 = p_helper;
      return(Rcpp::DataFrame::create(Rcpp::Named("ind") = ind, Rcpp::Named("prob1") = omega1, Rcpp::Named("prob2") = omega2, Rcpp::Named("theta21") = theta21, Rcpp::Named("theta22") = theta22,
                                     Rcpp::Named("loglik") = loglik[k-1], Rcpp::Named("iteration") = k,
                                     Rcpp::Named("censored_rate") = 1-sum(d)/d.size(), Rcpp::Named("message") = "convergent"));
    }else{
      double theta21 = -delta/std::log(1 - p1);
      double theta22 = -delta/std::log(1 - p2);
      return(Rcpp::DataFrame::create(Rcpp::Named("ind") = ind, Rcpp::Named("prob1") = omega1, Rcpp::Named("prob2") = omega2, Rcpp::Named("theta21") = theta21, Rcpp::Named("theta22") = theta22,
                                     Rcpp::Named("loglik") = loglik[k-1], Rcpp::Named("iteration") = k,
                                     Rcpp::Named("censored_rate") = 1-sum(d)/d.size(), Rcpp::Named("message") = "convergent"));
    }
  }
  else{
    return(Rcpp::DataFrame::create(Rcpp::Named("ind") = ind, Rcpp::Named("prob1") = R_NaN, Rcpp::Named("prob2") = R_NaN, Rcpp::Named("theta21") = R_NaN, Rcpp::Named("theta22") = R_NaN,
                                   Rcpp::Named("loglik") = R_NaN, Rcpp::Named("iteration") = k,
                                   Rcpp::Named("censored_rate") = 1-sum(d)/d.size(), Rcpp::Named("message") = "not convergent"));
  }
  
  
  
  
}