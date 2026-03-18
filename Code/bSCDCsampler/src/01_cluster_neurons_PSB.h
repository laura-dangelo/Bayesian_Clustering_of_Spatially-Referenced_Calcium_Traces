#ifndef PROBITSB_H
#define PROBITSB_H

#include "00_sample_mvnorm.h"
#include <truncnorm.h>
#include "00_SB.h"


Rcpp::List cluster_neurons(arma::mat AA, 
                           arma::mat logprobit_latentGP,
                           arma::vec cl1, int max_cl1,
                           const arma::mat inv_varcov_loc,
                           arma::mat latent_alpha,
                           double mu_alpha,
                           double sigma2_alpha) ;


                  
#endif

