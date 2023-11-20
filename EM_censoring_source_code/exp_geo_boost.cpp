#include "exp_geo_boost.h"


arma::vec dexp_boost(const arma::vec& x, const double& p) {
  boost::math::exponential dist(p);
  //boost::math::geometric_distribution<double> dist(p); // should work as well
  int n = x.size();
  arma::vec res(n);
  for (int i = 0; i < n; ++i) {
    res[i] = boost::math::pdf(dist, x[i]);
  }
  return res;
}


arma::vec pexp_boost(const arma::vec& x, const double& p) {
  boost::math::exponential dist(p);
  //boost::math::geometric_distribution<double> dist(p); // should work as well
  int n = x.size();
  arma::vec res(n);
  for (int i = 0; i < n; ++i) {
    res[i] = boost::math::cdf(dist, x[i]);
  }
  return res;
}

arma::vec dgeom_boost(const arma::vec& x, const double& p) {
  boost::math::geometric dist(p);
  //boost::math::geometric_distribution<double> dist(p); // should work as well
  int n = x.size();
  arma::vec res(n);
  for (int i = 0; i < n; ++i) {
    res[i] = boost::math::pdf(dist, x[i]);
  }
  return res;
}


arma::vec pgeom_boost(const arma::vec& x, const double& p) {
  boost::math::geometric dist(p);
  //boost::math::geometric_distribution<double> dist(p); // should work as well
  int n = x.size();
  arma::vec res(n);
  for (int i = 0; i < n; ++i) {
    res[i] = boost::math::cdf(dist, x[i]);
  }
  return res;
}

