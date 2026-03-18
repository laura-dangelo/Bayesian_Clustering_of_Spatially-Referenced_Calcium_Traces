#ifndef FUNS_H
#define FUNS_H

#include <RcppArmadillo.h>

arma::mat compute_Omega(int TT, double lGP_tau2_omega, double lGP_sigma2_omega) ;
arma::mat compute_Sigma(arma::mat loc, double theta) ;
arma::vec compute_cor(double t_new, int T, double lGP_tau2_omega, double lGP_sigma2_omega) ;

#endif

