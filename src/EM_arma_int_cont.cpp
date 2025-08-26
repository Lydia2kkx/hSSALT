// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]
#include <RcppArmadillo.h>
#include <boost/math/distributions/exponential.hpp>
#include <boost/math/distributions/geometric.hpp>

using namespace Rcpp;
using namespace arma;



arma::vec dexp_boost(arma::vec x, double p) {
  boost::math::exponential dist(p);
  //boost::math::geometric_distribution<double> dist(p); // should work as well
  int n = x.size();
  arma::vec res(n);
  for (int i = 0; i < n; ++i) {
    res[i] = boost::math::pdf(dist, x[i]);
  }
  return res;
}


arma::vec pexp_boost(arma::vec x, double p) {
  boost::math::exponential dist(p);
  //boost::math::geometric_distribution<double> dist(p); // should work as well
  int n = x.size();
  arma::vec res(n);
  for (int i = 0; i < n; ++i) {
    res[i] = boost::math::cdf(dist, x[i]);
  }
  return res;
}

// arma::vec mylog_cpp(arma::vec x) {
//   double small_offset = 1e-10;
//   int n = x.size();
//   arma::vec res(n);
//   for (int i = 0; i < n; ++i) {
//     if(x[i] <= 1e-10){
//       res[i] = log(x[i] + small_offset);
//     }else{
//       res[i] = log(x[i]);
//     }
//   }
//   return res;
// }

// arma::vec mylog_cpp(const arma::vec& x) {
//   const double small_offset = 1e-10;
//   return arma::log(arma::clamp(x, small_offset, arma::datum::inf) + small_offset);
// }


double mysum_cpp(arma::vec v ) {
  
  
  double s =0;
  
  int n = v.size();
  
  for (int i=0; i < n; i++) {
    // Guard against non-finite values
    if(is_finite(v[i])){
      s +=  v[i];
    }
    
  }
  return(s);
}


arma::vec dgeom_boost(arma::vec x, double p) {
  boost::math::geometric dist(p);
  //boost::math::geometric_distribution<double> dist(p); // should work as well
  int n = x.size();
  arma::vec res(n);
  for (int i = 0; i < n; ++i) {
    res[i] = boost::math::pdf(dist, x[i]);
  }
  return res;
}


arma::vec pgeom_boost(arma::vec x, double p) {
  boost::math::geometric dist(p);
  //boost::math::geometric_distribution<double> dist(p); // should work as well
  int n = x.size();
  arma::vec res(n);
  for (int i = 0; i < n; ++i) {
    res[i] = boost::math::cdf(dist, x[i]);
  }
  return res;
}



// [[Rcpp::export]]
DataFrame EM_algorithm_censored_arma(arma::vec data, double ind, arma::vec d, DataFrame parameter_starts, int N = 1000){
  DataFrame result;
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
  while ((abs(loglik[k]-loglik[k-1]) >= 1e-8) && (k <= (N-1))) {
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
  
  
  
  if(k <= N){
    double p_helper = pi1;
    if(rate1 <= rate2){
      double theta21 = 1/rate2;
      double theta22 = 1/rate1;
      double pi1 = pi2;
      double pi2 = p_helper;
      return(DataFrame::create(Named("ind") = ind, Named("pi1") = pi1, Named("pi2") = pi2, Named("theta21") = theta21, Named("theta22") = theta22,
                                     Named("loglik") = loglik[k-1], Named("iteration") = k,
                                     Named("censored_rate") = 1-sum(d)/d.size(), Named("message") = "convergent"));
    }else{
      double theta21 = 1/rate1;
      double theta22 = 1/rate2;
      return(DataFrame::create(Named("ind") = ind, Named("pi1") = pi1, Named("pi2") = pi2, Named("theta21") = theta21, Named("theta22") = theta22,
                                     Named("loglik") = loglik[k-1], Named("iteration") = k,
                                     Named("censored_rate") = 1-sum(d)/d.size(), Named("message") = "convergent"));
    }
  }
  else{
    return(DataFrame::create(Named("ind") = ind, Named("pi1") = R_NaN, Named("pi2") = R_NaN, Named("theta21") = R_NaN, Named("theta22") = R_NaN,
                                   Named("loglik") = R_NaN, Named("iteration") = k-1,
                                   Named("censored_rate") = 1-sum(d)/d.size(), Named("message") = "not convergent"));
  }
  
}


// [[Rcpp::export]]
DataFrame EM_algorithm_interval_arma(arma::vec data, double delta, double ind, arma::vec d, DataFrame parameter_starts, double q2, int N = 1000){
  arma::vec firstcol = parameter_starts[0];
  arma::vec secondcol = parameter_starts[1];
  arma::vec thirdcol = parameter_starts[2];
  
  double omega1 = firstcol[ind-1];
  
  double omega2 =  1 - omega1;
  double p1 = 1 - exp(-delta/secondcol[(ind-1)]);
  double p2 = 1 - exp(-delta/thirdcol[(ind-1)]);
  
  DataFrame result;
  bool Flag = false;
  
  // std::cout << "outside\n";
  // 
  // std::cout << "omega1 " << omega1 << "\n";
  // std::cout << "omega2 " << omega2<< "\n";
  // 
  // std::cout << "p1 " << p1 << "\n";
  // std::cout << "p2 " << p2 << "\n";
  
  
  
  arma::vec loglik(N);
  loglik[0] = 0;
  
  loglik[1]=mysum_cpp(omega1 * (log(omega1)+d%log(dgeom_boost(data, p1))+(1-d)%log(1-pgeom_boost(data, p1))))+ 
    mysum_cpp(omega2 * (log(omega2)+log(dgeom_boost(data, p2))+(1-d)%log(1-pgeom_boost(data, p2))));          
  
  int k =1;
  
  // loop
  while ((abs(loglik[k]-loglik[k-1]) >= 1e-8) && (k <= (N-1))) {
    // std::cout << "inside\n";
    // 
    // std::cout << "omega1 " << omega1 << "\n";
    // std::cout << "omega2 " << omega2<< "\n";
    // 
    // std::cout << "p1 " << p1 << "\n";
    // std::cout << "p2 " << p2 << "\n";
    
    // E step
    arma::vec sum_of_comps1 = omega1 * dgeom_boost(data, p1)+omega2 * dgeom_boost(data, p2); 
    arma::vec sum_of_comps2 = omega1 * (1-pgeom_boost(data, p1))+omega2 * (1-pgeom_boost(data, p2));  
    arma::vec tau1 = arma::pow((omega1 * dgeom_boost(data, p1)/sum_of_comps1),d)%arma::pow((omega1 * (1-pgeom_boost(data, p1))/sum_of_comps2),(1-d)); 
    arma::vec tau2 = arma::pow((omega2 * dgeom_boost(data, p2)/sum_of_comps1),d)%arma::pow((omega2 * (1-pgeom_boost(data, p2))/sum_of_comps2),(1-d)); 
    
    // std::cout << "sum_of_comps1 " << sum_of_comps1[0] << "\n";
    // std::cout << "sum_of_comps2 " << sum_of_comps2[0] << "\n";
    // 
    // std::cout << "tau1 " << tau1[0] << "\n";
    // std::cout << "tau2 " << tau2[0] << "\n";
    
    // M step
    omega1 = mysum_cpp(tau1)/mysum_cpp(tau1 + tau2); 
    omega2 = mysum_cpp(tau2)/mysum_cpp(tau1 + tau2); 
    
    p1 = mysum_cpp(tau1%d)/mysum_cpp(tau1%(d%(data+1) + (1-d) * q2));
    p2 = mysum_cpp(tau2%d)/mysum_cpp(tau2%(d%(data+1) + (1-d) * q2));
    if (std::abs(p1 - 1) < 1e-128 || std::abs(p2 - 1) < 1e-128) {
      // std::cout << "Warning: p1 reached 1. Exiting EM algorithm early." << std::endl;
      double p_helper = omega1;
      if(p1 <= p2){
        double theta21 = 0;
        double theta22 = -delta/std::log(1 - p1);
        double omega1 = omega2;
        double omega2 = p_helper;
        result = DataFrame::create(Named("ind") = ind, Named("prob1") = omega1, Named("prob2") = omega2, Named("theta21") = theta21, Named("theta22") = theta22,
                                         Named("loglik") = loglik[k-1], Named("iteration") = k,
                                         Named("censored_rate") = 1-sum(d)/d.size(), Named("message") = "p1 reached 1");
      }else{
        double theta21 = 0;
        double theta22 = -delta/std::log(1 - p2);
        result = DataFrame::create(Named("ind") = ind, Named("prob1") = omega1, Named("prob2") = omega2, Named("theta21") = theta21, Named("theta22") = theta22,
                                         Named("loglik") = loglik[k-1], Named("iteration") = k,
                                         Named("censored_rate") = 1-sum(d)/d.size(), Named("message") = "p1 reached 1");
      }
      Flag = true;
      break;
    }
    
    // std::cout << "omega1 " << omega1 << "\n";
    // std::cout << "omega2 " << omega2<< "\n";
    // 
    // std::cout << "p1 " << p1 << "\n";
    // std::cout << "p2 " << p2 << "\n";
    
    loglik[k+1]=mysum_cpp(tau1%(log(omega1)+d%arma::log(dgeom_boost(data, p1)) + (1-d)%arma::log(1-pgeom_boost(data, p1))))+ 
      mysum_cpp(tau2%(log(omega2)+d%arma::log(dgeom_boost(data, p2))+ (1-d)%arma::log(1-pgeom_boost(data, p2))));           
    
    // std::cout << "loglik k+1 " << loglik[k+1] << "\n";
    
    k++;
    
    // std::cout << "k+1 = " << k << "\n";
  }
  if (!Flag) {
    if(k <= N){
      double p_helper = omega1;
      if(p1 <= p2){
        double theta21 = -delta/std::log(1 - p2);
        double theta22 = -delta/std::log(1 - p1);
        double omega1 = omega2;
        double omega2 = p_helper;
        result = DataFrame::create(Named("ind") = ind, Named("prob1") = omega1, Named("prob2") = omega2, Named("theta21") = theta21, Named("theta22") = theta22,
                                         Named("loglik") = loglik[k-1], Named("iteration") = k,
                                         Named("censored_rate") = 1-sum(d)/d.size(), Named("message") = "convergent");
      }else{
        double theta21 = -delta/std::log(1 - p1);
        double theta22 = -delta/std::log(1 - p2);
        result = DataFrame::create(Named("ind") = ind, Named("prob1") = omega1, Named("prob2") = omega2, Named("theta21") = theta21, Named("theta22") = theta22,
                                         Named("loglik") = loglik[k-1], Named("iteration") = k,
                                         Named("censored_rate") = 1-sum(d)/d.size(), Named("message") = "convergent");
      }
    }else{
      result = DataFrame::create(Named("ind") = ind, Named("prob1") = R_NaN, Named("prob2") = R_NaN, Named("theta21") = R_NaN, Named("theta22") = R_NaN,
                                       Named("loglik") = R_NaN, Named("iteration") = k-1,
                                       Named("censored_rate") = 1-sum(d)/d.size(), Named("message") = "not convergent");
    }
  }
  return result;
  
}
