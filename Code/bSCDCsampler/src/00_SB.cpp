#include "00_SB.h"
#include <RcppArmadilloExtensions/sample.h>

/*
 * compute stick-breaking starting from beta r.v.
 */ 
arma::vec stick_breaking(arma::vec beta_var)
{
  int len = beta_var.n_elem ;
  arma::vec out(len) ;
  
  out(0) = beta_var(0) ;
  for(int k = 1; k < len; k++)
  {
    out(k) = beta_var(k) * arma::prod( 1.0 - beta_var.head(k) ) ;
  }
  return(out) ;
}


arma::vec log_stick_breaking(arma::vec beta_var)
{
  int len = beta_var.n_elem ;
  arma::vec one_m_beta = arma::ones(len) - beta_var ;
  arma::vec out(len) ;
  
  out(0) = log(beta_var(0)) ;
  for(int k = 1; k < len; k++)
  {
    out(k) = log(beta_var(k)) + arma::accu( log(one_m_beta.head(k)) ) ;
  }
  return(out) ;
}


int sample_cl(arma::vec cluster_id, arma::vec prob) {
  int out ; 
  out = Rcpp::RcppArmadillo::sample(cluster_id, 1, false, prob)[0] ;
  return out ;
}