#include "misc.h"


// using namespace Rcpp;
// using namespace arma;


Rcpp::IntegerVector table_factor_cpp(const Rcpp::IntegerVector& v){
  Rcpp::CharacterVector ch = v.attr("levels");
  Rcpp::IntegerVector v2(ch.size());
  Rcpp::IntegerVector v3(ch.size());
  for (int i = 0; i != ch.size(); ++i) {
    v2[i]= i+1;
    v3[i]= std::count(v.begin(),v.end(),v2.at(i));
  }
  v3.attr("names") = ch;
  return v3;
}

// Sort from this internet source
// https://gallery.rcpp.org/articles/sorting/


Rcpp::NumericVector stl_sort(const Rcpp::NumericVector& x){
  Rcpp::NumericVector y = clone(x);
  std::sort(y.begin(), y.end());
  return y;
}



double mysum_cpp(const arma::vec& v ) {
  
  double s =0;
  int n = v.size();
  
  for (int i=0; i < n; i++) {
    // Guard against non-finite values
    if(arma::is_finite(v[i])){
      s +=  v[i];
    }
  }
  return(s);
}



Rcpp::NumericVector rexpC(const double& n  , const double& par) {
  return Rcpp::rexp( n, par); 
}


