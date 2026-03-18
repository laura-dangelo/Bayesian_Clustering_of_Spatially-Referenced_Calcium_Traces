#ifndef STATESPACE_H
#define STATESPACE_H

#include "00_sample_mvnorm.h"
#include <truncnorm.h>
                
arma::mat sample_latent_signal(arma::mat AA, 
                              arma::mat latent_signal, 
                              arma::vec cl_s,
                              int max_cl,
                              double lGP_tau2_omega, 
                              double lGP_sigma2_omega, 
                              double lGP_error_variance, 
                              arma::mat Omega,
                              arma::mat invOmega,
                              arma::mat kk) ;


#endif

