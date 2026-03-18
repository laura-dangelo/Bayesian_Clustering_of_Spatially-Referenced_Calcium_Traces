#ifndef MHGAMMA_H
#define MHGAMMA_H

#include <RcppArmadillo.h>

double logprior_gamma(double gamma, double hyp_gamma1, double hyp_gamma2) ;

double logpost_gamma(arma::mat Ca, arma::mat AA, 
                     double gamma, double tau2, 
                     double hyp_gamma1, double hyp_gamma2) ;

                  
#endif

