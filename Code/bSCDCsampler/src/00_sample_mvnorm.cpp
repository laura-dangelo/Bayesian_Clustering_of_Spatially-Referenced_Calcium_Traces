#include "00_sample_mvnorm.h"


/*
 * sample multivariate normal distribution
 */
arma::mat rmvnorm(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}

