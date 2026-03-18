#ifndef SB_H
#define SB_H

#include <RcppArmadillo.h>


arma::vec stick_breaking(arma::vec beta_var) ;
arma::vec log_stick_breaking(arma::vec beta_var) ;
int sample_cl(arma::vec cluster_id, arma::vec prob) ;

#endif

