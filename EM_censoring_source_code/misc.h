#ifndef MISC_H
#define MISC_H

#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]



// [[Rcpp::export]] 
Rcpp::IntegerVector table_factor_cpp(const Rcpp::IntegerVector& v);

// [[Rcpp::export]] 
Rcpp::NumericVector stl_sort(const Rcpp::NumericVector& x);

// [[Rcpp::export]] 
double mysum_cpp(const arma::vec& v );

// [[Rcpp::export]] 
Rcpp::NumericVector rexpC(const double& n  , const double& par);


#endif