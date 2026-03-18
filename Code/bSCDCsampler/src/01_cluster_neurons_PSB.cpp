#include "01_cluster_neurons_PSB.h"


/*
 * cluster allocation using location-dependent probit-stick-breaking process
 */
Rcpp::List cluster_neurons(arma::mat AA, 
                           arma::mat logprobit_latentGP,
                           arma::vec cl_s, 
                           int max_cl,
                           const arma::mat inv_varcov_loc,
                           arma::mat latent_alpha,
                           double mu_alpha,
                           double sigma2_alpha)
{
  int n = AA.n_cols ;
  int T = AA.n_rows ;
  bool error = false ;

  arma::vec w1(max_cl) ;
  arma::vec prob1(max_cl) ;
  arma::mat latent_z(n, max_cl) ;
  
  arma::vec tmp_lik_prob(max_cl) ;
  arma::vec tmp_probit_alpha(max_cl) ;
  
  // reconstruct the binary signal for each (i,t)
  // matrix T x n, binary_signal(t,i) = I( a(t,i) > 0 )
  arma::mat binary_signal(T, n) ;
  binary_signal.fill(0.0) ;
  for(int i = 0; i < n; i++) {
    for(int t = 0; t < T; t++) {
      if(AA(t,i) > 0.0) { binary_signal(t,i) = 1.0 ; }
    }
  }
  
  
  arma::vec tmp(T) ;
  
  /*
   THIS IS STEP (a) IN THE ALGORITHM
   Compute stick breaking based on the data augmentation for PSB (1)
   compute the posterior cluster probabilities (based on z, alpha, and likelihood) (2)
   sample cluster allocation (3)
   for each i = 1,...,n 
   */

  // for loop to compute the probabilities and allocate one neuron at a time
  for(int i = 0; i < n; i++) {
    
    for(int k = 0; k < max_cl; k++) { 
      // probability of binary signal s_it given the k-th atom (latent GP \tilde{s}*_k)
      tmp = binary_signal.col(i) % logprobit_latentGP.col(k) + 
        ( arma::ones(T) - binary_signal.col(i) ) % log( arma::ones(T) + arma::vec(T, arma::fill::value(1e-10)) - exp(logprobit_latentGP.col(k)) ) ;
      
      tmp_lik_prob(k) = arma::accu( tmp ) ;
      
      // prior probability of the PSBP
      // latent_alpha is an (n x max_cl) matrix
      tmp_probit_alpha(k) = R::pnorm( latent_alpha(i,k), 0.0, 1.0, true, false ) ;
      
    }
    
    // PSBP weights from the probit-transformed latent_alpha
    w1 = log(stick_breaking(tmp_probit_alpha)) ;
    
    // combine the two probabilities
    prob1 = tmp_lik_prob + w1 ; 
    prob1 =  prob1 - max(prob1) ;
    prob1 = exp(prob1) ;
    
    if(!arma::is_finite(prob1)) {
      prob1.fill(0.0) ;
      prob1(0) = 1.0 ;
      error = true ;  // Rcpp::Rcout << "error sampling loc_dep_kernelSB" ;
    }
    
    arma::vec cluster_id = arma::linspace(0, max_cl-1, max_cl) ;
    cl_s(i) = sample_cl(cluster_id, prob1) ;
  }
  
  
  
  
  /*
   * THIS IS STEP (b)
   */
  /*
   * Step (b1)
   * Sample the latent truncated Gaussian random variables Z (data augmentation PSBP)
   */
  latent_z.zeros() ;
  for(int i = 0; i < n; i++)
  {
    // if k < (cluster allocation of neuron i)  => generate from a truncated normal between (-infty, 0)
    if(cl_s(i) > 0) {
      for(int k = 0; k < cl_s(i); k++) {
        latent_z(i, k) = r_truncnorm(latent_alpha(i, k), 1, -999, 0) ; // truncated on (-infty, 0)
      }
    }
    // if k = (cluster allocation of neuron i)  => generate from a truncated normal between (0, +infty)
    latent_z(i, cl_s(i)) = r_truncnorm(latent_alpha(i, cl_s(i)), 1, 0, 999) ; // truncated on (0, +infty)
  }
  
  /*
   * Step (b2)
   * Sample the latent multivariate Gaussian r.v. alpha (with loc-dep covariance matrix) 
   */
  arma::mat par2 = inv_sympd(inv_varcov_loc +  (1.0 / sigma2_alpha) * arma::eye(n,n)) ; 
  for(int k = 0; k < max_cl; k++) {
    latent_alpha.col(k) = rmvnorm(1, par2 * (mu_alpha * inv_varcov_loc * arma::ones(n) + (1.0 / sigma2_alpha) * latent_z.col(k)), par2).t() ;
  }

 
  
  return Rcpp::List::create(Rcpp::Named("latent_alpha") = latent_alpha,
                            Rcpp::Named("cluster_signal") = cl_s,
                            Rcpp::Named("error") = error);
}

