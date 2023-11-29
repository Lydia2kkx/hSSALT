#ifndef EM_ALGORITHM_CONT_H
#define EM_ALGORITHM_CONT_H

#include <RcppArmadillo.h>
#include <Rcpp.h>
// #include <boost/math/distributions/exponential.hpp>
// #include <boost/math/distributions/geometric.hpp>

#include "misc.h"
#include "exp_geo_boost.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]

// [[Rcpp::export]]
Rcpp::List EM_algorithm_censored_arma(const arma::vec& data, const double& ind, const arma::vec& d, const Rcpp::DataFrame& parameter_starts, const int& N, const double& tol);

#endif
