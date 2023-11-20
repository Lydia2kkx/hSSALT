# include "EM_algorithm_cont.h"




Rcpp::DataFrame EM_algorithm_censored_arma(const arma::vec& data, const double& ind, const arma::vec& d, const Rcpp::DataFrame& parameter_starts, const int& N, const double& tol){
  arma::vec firstcol = parameter_starts[0];
  arma::vec secondcol = parameter_starts[1];
  arma::vec thirdcol = parameter_starts[2];
  
  double pi1 = firstcol[ind-1];
  
  double pi2 =  1 - pi1;
  double rate1 = 1/secondcol[(ind-1)];
  double rate2 = 1/thirdcol[(ind-1)];
  
  // std::cout << "outside\n";
  // 
  // std::cout << "pi1 " << pi1 << "\n";
  // std::cout << "pi2 " << pi2<< "\n";
  // 
  // std::cout << "rate1 " << rate1 << "\n";
  // std::cout << "rate2 " << rate2 << "\n";
  
  
  
  arma::vec loglik(N);
  loglik[0] = 0;
  
  loglik[1]=mysum_cpp(pi1 * (log(pi1)+d%log(dexp_boost(data, rate1))+(1-d)%log(1-pexp_boost(data, rate1))))+ 
    mysum_cpp(pi2 * (log(pi2)+log(dexp_boost(data, rate2))+(1-d)%log(1-pexp_boost(data, rate2)))); 
  
  // std::cout << "loglik1 " << loglik[1] << "\n";
  
  int k =1;
  
  // loop
  while ((abs(loglik[k]-loglik[k-1]) >= tol) && (k <= (N-1))) {
    // std::cout << "inside\n";
    // 
    // std::cout << "pi1 " << pi1 << "\n";
    // std::cout << "pi2 " << pi1 << "\n";
    // 
    // std::cout << "rate1 " << rate1 << "\n";
    // std::cout << "rate2 " << rate2 << "\n";
    
    // E step
    arma::vec sum_of_comps1 = pi1 * dexp_boost(data, rate1)+pi2 * dexp_boost(data, rate2); 
    arma::vec sum_of_comps2 = pi1 * (1-pexp_boost(data, rate1))+pi2 * (1-pexp_boost(data, rate2));  
    arma::vec tau1 = arma::pow((pi1 * dexp_boost(data, rate1)/sum_of_comps1),d)%arma::pow((pi1 * (1-pexp_boost(data, rate1))/sum_of_comps2),(1-d)); 
    arma::vec tau2 = arma::pow((pi2 * dexp_boost(data, rate2)/sum_of_comps1),d)%arma::pow((pi2 * (1-pexp_boost(data, rate2))/sum_of_comps2),(1-d)); 
    
    // std::cout << "sum_of_comps1 " << sum_of_comps1[0] << "\n";
    // std::cout << "sum_of_comps2 " << sum_of_comps2[0] << "\n";
    //  
    // std::cout << "tau1 " << tau1[0] << "\n";
    // std::cout << "tau2 " << tau2[0] << "\n";
    
    // M step
    pi1 = mysum_cpp(tau1)/data.size(); 
    pi2 = mysum_cpp(tau2)/data.size(); 
    
    rate1 = mysum_cpp(tau1%d)/mysum_cpp(tau1%(data));
    rate2 = mysum_cpp(tau2%d)/mysum_cpp(tau2%(data));
    
    // std::cout << "pi1 " << pi1 << "\n";
    // std::cout << "pi2 " << pi2<< "\n";
    // 
    // std::cout << "rate1 " << rate1 << "\n";
    // std::cout << "rate2 " << rate2 << "\n";
    
    loglik[k+1]=mysum_cpp(tau1%(log(pi1)+d%arma::log(dexp_boost(data, rate1)) + (1-d)%arma::log(1-pexp_boost(data, rate1))))+ 
      mysum_cpp(tau2%(log(pi2)+d%arma::log(dexp_boost(data, rate2))+ (1-d)%arma::log(1-pexp_boost(data, rate2))));           
    
    
    // std::cout << "loglik k+1 " << loglik[k+1] << "\n";
    
    k++;
    
    // std::cout << "k+1 = " << k << "\n";
  }
  
  
  
  if(k <= (N-1)){
    double p_helper = pi1;
    if(rate1 <= rate2){
      double theta21 = 1/rate2;
      double theta22 = 1/rate1;
      double pi1 = pi2;
      double pi2 = p_helper;
      return(Rcpp::DataFrame::create(Rcpp::Named("ind") = ind, Rcpp::Named("pi1") = pi1, Rcpp::Named("pi2") = pi2, Rcpp::Named("theta21") = theta21, Rcpp::Named("theta22") = theta22,
                                     Rcpp::Named("loglik") = loglik[k-1], Rcpp::Named("iteration") = k,
                                     Rcpp::Named("censored_rate") = 1-sum(d)/d.size(), Rcpp::Named("message") = "convergent"));
    }else{
      double theta21 = 1/rate1;
      double theta22 = 1/rate2;
      return(Rcpp::DataFrame::create(Rcpp::Named("ind") = ind, Rcpp::Named("pi1") = pi1, Rcpp::Named("pi2") = pi2, Rcpp::Named("theta21") = theta21, Rcpp::Named("theta22") = theta22,
                                     Rcpp::Named("loglik") = loglik[k-1], Rcpp::Named("iteration") = k,
                                     Rcpp::Named("censored_rate") = 1-sum(d)/d.size(), Rcpp::Named("message") = "convergent"));
    }
  }
  else{
    return(Rcpp::DataFrame::create(Rcpp::Named("ind") = ind, Rcpp::Named("pi1") = R_NaN, Rcpp::Named("pi2") = R_NaN, Rcpp::Named("theta21") = R_NaN, Rcpp::Named("theta22") = R_NaN,
                                   Rcpp::Named("loglik") = R_NaN, Rcpp::Named("iteration") = k,
                                   Rcpp::Named("censored_rate") = 1-sum(d)/d.size(), Rcpp::Named("message") = "not convergent"));
  }
  
}
