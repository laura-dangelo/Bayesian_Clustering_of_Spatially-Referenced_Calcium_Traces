#ifndef AMPLITUDES_H
#define AMPLITUDES_H

#include "00_SB.h"


Rcpp::List cluster_amplitudes(arma::mat calcium,
                              double gamma,
                              arma::mat latent_signal, 
                              arma::vec cl_s,
                              double tau2,
                              arma::vec SB_weights,
                              arma::vec unique_ampl,
                              arma::mat cl_a,
                              int max_cl,
                              double alpha,
                              double hyp_a1, double hyp_a2,
                              arma::vec eps_a_v, 
                              double step_AMH_a,
                              int iter,
                              arma::mat acc_vec_MHa,
                              double a_low
                              ) ;
                  
#endif
