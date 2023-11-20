#ifndef EXP_GEO_BOOST_H
#define EXP_GEO_BOOST_H

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <boost/math/distributions/exponential.hpp>
#include <boost/math/distributions/geometric.hpp>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]

// [[Rcpp::export]] 
arma::vec dexp_boost(const arma::vec& x, const double& p);

// [[Rcpp::export]] 
arma::vec pexp_boost(const arma::vec& x, const double& p);

// [[Rcpp::export]] 
arma::vec dgeom_boost(const arma::vec& x, const double& p);

// [[Rcpp::export]] 
arma::vec pgeom_boost(const arma::vec& x, const double& p);

# endif 

