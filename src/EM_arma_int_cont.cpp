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
List EM_algorithm_censored_arma(arma::vec data, double ind, arma::vec d, DataFrame parameter_starts, int N = 1000, double tol = 1e-8){
  DataFrame results;
  arma::vec firstcol = parameter_starts[0];
  arma::vec secondcol = parameter_starts[1];
  arma::vec thirdcol = parameter_starts[2];
  
  double p1 = firstcol[ind-1];
  
  double p2 =  1 - p1;
  double rate1 = 1/secondcol[(ind-1)];
  double rate2 = 1/thirdcol[(ind-1)];
  
  arma::vec loglik(N);
  loglik[0] = 0;
  
  loglik[1]=mysum_cpp(p1 * (log(p1)+d%log(dexp_boost(data, rate1))+(1-d)%log(1-pexp_boost(data, rate1))))+ 
    mysum_cpp(p2 * (log(p2)+log(dexp_boost(data, rate2))+(1-d)%log(1-pexp_boost(data, rate2)))); 
  
  int k = 1;
  
  arma::vec tau1;
  arma::vec tau2; 
  
  // loop
  while ((abs(loglik[k]-loglik[k-1]) >= tol) && (k <= (N-1))) {
    
    // E step
    arma::vec sum_of_comps1 = p1 * dexp_boost(data, rate1)+p2 * dexp_boost(data, rate2); 
    arma::vec sum_of_comps2 = p1 * (1-pexp_boost(data, rate1))+p2 * (1-pexp_boost(data, rate2));  
    
    tau1 = arma::pow((p1 * dexp_boost(data, rate1)/sum_of_comps1),d)%arma::pow((p1 * (1-pexp_boost(data, rate1))/sum_of_comps2),(1-d)); 
    tau2 = arma::pow((p2 * dexp_boost(data, rate2)/sum_of_comps1),d)%arma::pow((p2 * (1-pexp_boost(data, rate2))/sum_of_comps2),(1-d)); 
    
    // M step
    p1 = mysum_cpp(tau1)/data.size(); 
    p2 = mysum_cpp(tau2)/data.size(); 
    
    rate1 = mysum_cpp(tau1%d)/mysum_cpp(tau1%(data));
    rate2 = mysum_cpp(tau2%d)/mysum_cpp(tau2%(data));
    
    loglik[k+1]=mysum_cpp(tau1%(log(p1)+d%arma::log(dexp_boost(data, rate1)) + (1-d)%arma::log(1-pexp_boost(data, rate1))))+ 
      mysum_cpp(tau2%(log(p2)+d%arma::log(dexp_boost(data, rate2))+ (1-d)%arma::log(1-pexp_boost(data, rate2))));           
    
    k++;
  
  }
  
  double p_helper = p1;
  double theta21, theta22;
  if(rate1 <= rate2){
    theta21 = 1/rate2;
    theta22 = 1/rate1;
    p1 = p2;
    p2 = p_helper;
  }else{
    theta21 = 1/rate1;
    theta22 = 1/rate2;
  }
  
  std::string message;
  arma::mat posterior;
  // Avner: Changed from k<=N because of differences in lists in C++ and R. Previously never went into the 'else' part
  if(k <= (N-1)){
    message = "convergent";
    posterior.set_size(data.n_elem, 2);
    posterior.col(0) = tau1;
    posterior.col(1) = tau2;
  } else{
    message = "not convergent";
    posterior.set_size(0,0);
  }
  
  results = DataFrame::create(Named("p1") = p1, Named("p2") = p2, Named("theta21") = theta21, Named("theta22") = theta22,
                                    Named("loglik") = loglik[k], Named("iteration") = k, Named("message") = message);
  
  return List::create(
    _["results"] = results,
    _["posterior"] = posterior
  );
  
}


// [[Rcpp::export]]
List EM_algorithm_interval_arma(arma::vec data, double delta, double ind, arma::vec d, DataFrame parameter_starts, double q2, int N = 1000, double tol = 1e-8){
  arma::vec firstcol = parameter_starts[0];
  arma::vec secondcol = parameter_starts[1];
  arma::vec thirdcol = parameter_starts[2];
  
  double omega1 = firstcol[ind-1];
  
  double omega2 =  1 - omega1;
  double p1 = 1 - exp(-delta/secondcol[(ind-1)]);
  double p2 = 1 - exp(-delta/thirdcol[(ind-1)]);
  
  DataFrame results;
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
  
  int k = 1;
  
  arma::vec tau1;
  arma::vec tau2; 
  
  // loop
  while ((abs(loglik[k]-loglik[k-1]) >= tol) && (k <= (N-1))) {
    
    // E step
    arma::vec sum_of_comps1 = omega1 * dgeom_boost(data, p1)+omega2 * dgeom_boost(data, p2); 
    arma::vec sum_of_comps2 = omega1 * (1-pgeom_boost(data, p1))+omega2 * (1-pgeom_boost(data, p2));  
    tau1 = arma::pow((omega1 * dgeom_boost(data, p1)/sum_of_comps1),d)%arma::pow((omega1 * (1-pgeom_boost(data, p1))/sum_of_comps2),(1-d)); 
    tau2 = arma::pow((omega2 * dgeom_boost(data, p2)/sum_of_comps1),d)%arma::pow((omega2 * (1-pgeom_boost(data, p2))/sum_of_comps2),(1-d)); 
    
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
        results = DataFrame::create(Named("p1") = omega1, Named("p2") = omega2, Named("theta21") = theta21, Named("theta22") = theta22,
                                         Named("loglik") = loglik[k-1], Named("iteration") = k, Named("message") = "p1 reached 1");
        //Avner: For now
        arma::mat posterior;
        posterior.set_size(0,0);
        return List::create(
          _["results"] = results,
          _["posterior"] = posterior
        );
      }else{
        double theta21 = 0;
        double theta22 = -delta/std::log(1 - p2);
        results = DataFrame::create(Named("p1") = omega1, Named("p2") = omega2, Named("theta21") = theta21, Named("theta22") = theta22,
                                         Named("loglik") = loglik[k-1], Named("iteration") = k, Named("message") = "p1 reached 1");
        //Avner: For now
        arma::mat posterior;
        posterior.set_size(0,0);
        return List::create(
          _["results"] = results,
          _["posterior"] = posterior
        );
      }
      Flag = true;
      break;
    }
    
    loglik[k+1]=mysum_cpp(tau1%(log(omega1)+d%arma::log(dgeom_boost(data, p1)) + (1-d)%arma::log(1-pgeom_boost(data, p1))))+ 
      mysum_cpp(tau2%(log(omega2)+d%arma::log(dgeom_boost(data, p2))+ (1-d)%arma::log(1-pgeom_boost(data, p2))));           
    
    k++;
  }
  
  
  if (!Flag) {
    double p_helper = omega1;
    double theta21, theta22;
    if(p1 <= p2){
      theta21 = -delta/std::log(1 - p2);
      theta22 = -delta/std::log(1 - p1);
      omega1 = omega2;
      omega2 = p_helper;
    }else{
      theta21 = -delta/std::log(1 - p1);
      theta22 = -delta/std::log(1 - p2);
    }
    
    std::string message;
    arma::mat posterior;
    if(k <= (N-1)){
      message = "convergent";
      posterior.set_size(data.n_elem, 2);
      posterior.col(0) = tau1;
      posterior.col(1) = tau2;
    } else{
      message = "not convergent";
      posterior.set_size(0,0);
    }
    
    results = DataFrame::create(Named("p1") = omega1, Named("p2") = omega2, Named("theta21") = theta21, Named("theta22") = theta22,
                                      Named("loglik") = loglik[k], Named("iteration") = k, Named("message") = message);
    
    return List::create(
      _["results"] = results,
      _["posterior"] = posterior
    );
    
  }
  
}
