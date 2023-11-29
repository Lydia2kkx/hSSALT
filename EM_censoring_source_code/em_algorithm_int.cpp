#include "em_algorithm_int.h"



Rcpp::List EM_algorithm_interval_arma(const arma::vec& data, const double& delta, const double& ind, const arma::vec& d, const Rcpp::DataFrame& parameter_starts, const double& q2, const int& N, const double& tol){
  arma::vec firstcol = parameter_starts[0];
  arma::vec secondcol = parameter_starts[1];
  arma::vec thirdcol = parameter_starts[2];

  double omega1 = firstcol[ind-1];

  double omega2 =  1 - omega1;
  double p1 = 1 - exp(-delta/secondcol[(ind-1)]);
  double p2 = 1 - exp(-delta/thirdcol[(ind-1)]);

  arma::vec tau1;
  arma::vec tau2;


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
    tau1 = arma::pow((omega1 * dgeom_boost(data, p1)/sum_of_comps1),d)%arma::pow((omega1 * (1-pgeom_boost(data, p1))/sum_of_comps2),(1-d));
    tau2 = arma::pow((omega2 * dgeom_boost(data, p2)/sum_of_comps1),d)%arma::pow((omega2 * (1-pgeom_boost(data, p2))/sum_of_comps2),(1-d));



    // M step
    omega1 = mysum_cpp(tau1)/mysum_cpp(tau1 + tau2);
    omega2 = mysum_cpp(tau2)/mysum_cpp(tau1 + tau2);

    p1 = mysum_cpp(tau1%d)/mysum_cpp(tau1%(d%(data+1) + (1-d) * q2));
    p2 = mysum_cpp(tau2%d)/mysum_cpp(tau2%(d%(data+1) + (1-d) * q2));



    loglik[k+1]=mysum_cpp(tau1%(log(omega1)+d%arma::log(dgeom_boost(data, p1)) + (1-d)%arma::log(1-pgeom_boost(data, p1))))+
      mysum_cpp(tau2%(log(omega2)+d%arma::log(dgeom_boost(data, p2))+ (1-d)%arma::log(1-pgeom_boost(data, p2))));



    k++;


  }

  double theta21;
  double theta22;
  // Rcpp::DataFrame posterior;
  Rcpp::NumericMatrix posterior;

  if(k <= (N-1)){
    double p_helper = omega1;
    if(p1 <= p2){
      theta21 = -delta/std::log(1 - p2);
      theta22 = -delta/std::log(1 - p1);
      omega1 = omega2;
      omega2 = p_helper;

    }else{
      theta21 = -delta/std::log(1 - p1);
      theta22 = -delta/std::log(1 - p2);

    }
    // posterior = Rcpp::DataFrame::create(Rcpp::Named("tau1")=tau1, Rcpp::Named("tau2")=tau2);
    posterior = Rcpp::wrap(arma::join_horiz(tau1,tau2));
    colnames(posterior) = Rcpp::CharacterVector::create("tau1", "tau2");

    return(Rcpp::List::create(Rcpp::Named("results")=Rcpp::DataFrame::create(Rcpp::Named("ind") = ind, Rcpp::Named("p1") = omega1, Rcpp::Named("p2") = omega2, Rcpp::Named("theta21") = theta21, Rcpp::Named("theta22") = theta22,
                                               Rcpp::Named("loglik") = loglik[k-1], Rcpp::Named("iteration") = k,
                                               Rcpp::Named("message") = "convergent"), Rcpp::Named("posterior") = posterior));
  }
  else{
    return(Rcpp::List::create(Rcpp::Named("results")=Rcpp::DataFrame::create(Rcpp::Named("ind") = ind, Rcpp::Named("p1") = R_NaN, Rcpp::Named("p2") = R_NaN, Rcpp::Named("theta21") = R_NaN, Rcpp::Named("theta22") = R_NaN,
                                   Rcpp::Named("loglik") = R_NaN, Rcpp::Named("iteration") = k,
                                   Rcpp::Named("message") = "not convergent"),Rcpp::Named("posterior") = R_NaN));
  }




}
