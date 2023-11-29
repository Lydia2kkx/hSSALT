# include "EM_algorithm_cont.h"




Rcpp::List EM_algorithm_censored_arma(const arma::vec& data, const double& ind, const arma::vec& d, const Rcpp::DataFrame& parameter_starts, const int& N, const double& tol){
  arma::vec firstcol = parameter_starts[0];
  arma::vec secondcol = parameter_starts[1];
  arma::vec thirdcol = parameter_starts[2];

  double p1 = firstcol[ind-1];

  double p2 =  1 - p1;
  double rate1 = 1/secondcol[(ind-1)];
  double rate2 = 1/thirdcol[(ind-1)];

  arma::vec tau1;
  arma::vec tau2;


  arma::vec loglik(N);
  loglik[0] = 0;

  loglik[1]=mysum_cpp(p1 * (log(p1)+d%log(dexp_boost(data, rate1))+(1-d)%log(1-pexp_boost(data, rate1))))+
    mysum_cpp(p2 * (log(p2)+log(dexp_boost(data, rate2))+(1-d)%log(1-pexp_boost(data, rate2))));


  int k =1;

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

  double theta21;
  double theta22;
  // Rcpp::DataFrame posterior;
  // posterior = Rcpp::DataFrame::create(Rcpp::Named("tau1")=tau1, Rcpp::Named("tau2")=tau2);

  Rcpp::NumericMatrix posterior;


  if(k <= (N-1)){
    double p_helper = p1;
    if(rate1 <= rate2){
      theta21 = 1/rate2;
      theta22 = 1/rate1;
      p1 = p2;
      p2 = p_helper;


    }else{
      theta21 = 1/rate1;
      theta22 = 1/rate2;

    }
    posterior = Rcpp::wrap(arma::join_horiz(tau1,tau2));
    colnames(posterior) = Rcpp::CharacterVector::create("tau1", "tau2");
    return(Rcpp::List::create(Rcpp::Named("results")=Rcpp::DataFrame::create(Rcpp::Named("ind") = ind, Rcpp::Named("p1") = p1, Rcpp::Named("p2") = p2, Rcpp::Named("theta21") = theta21, Rcpp::Named("theta22") = theta22,
                                          Rcpp::Named("loglik") = loglik[k-1], Rcpp::Named("iteration") = k,
                                          Rcpp::Named("message") = "convergent"),Rcpp::Named("posterior") = posterior));
  }
  else{
    return(Rcpp::List::create(Rcpp::Named("results")=Rcpp::DataFrame::create(Rcpp::Named("ind") = ind, Rcpp::Named("p1") = R_NaN, Rcpp::Named("p2") = R_NaN, Rcpp::Named("theta21") = R_NaN, Rcpp::Named("theta22") = R_NaN,
                                   Rcpp::Named("loglik") = R_NaN, Rcpp::Named("iteration") = k,
                                   Rcpp::Named("message") = "not convergent"),Rcpp::Named("posterior") = R_NaN ));
  }

}
