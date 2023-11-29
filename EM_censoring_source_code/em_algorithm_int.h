#ifndef EM_ALGORITHM_INT_H
#define EM_ALGORITHM_INT_H

#include <RcppArmadillo.h>
#include <Rcpp.h>
// #include <boost/math/distributions/exponential.hpp>
// #include <boost/math/distributions/geometric.hpp>

#include "misc.h"
#include "exp_geo_boost.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]

// [[Rcpp::export]]
Rcpp::List EM_algorithm_interval_arma(const arma::vec& data, const double& delta, const double& ind, const arma::vec& d, const Rcpp::DataFrame& parameter_starts, const double& q2, const int& N, const double& tol);

#endif
