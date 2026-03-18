#ifndef NORMALRVFUNCTIONS_H
#define NORMALRVFUNCTIONS_H

#include <RcppArmadillo.h>

const double log2pi = std::log(2.0 * M_PI) ;

arma::mat rmvnorm(int n, arma::vec mu, arma::mat sigma) ;
                   
#endif

