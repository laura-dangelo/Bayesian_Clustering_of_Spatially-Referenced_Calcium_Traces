#include "01_sample_latent_signal.h"

/*
 * Atoms of the mixture on the latent signal: realizations of the GP
 */

arma::mat sample_latent_signal(arma::mat AA, 
                               arma::mat latent_signal, // GP realizations
                               arma::vec cl_s,
                               int max_cl,
                               double lGP_tau2_omega, 
                               double lGP_sigma2_omega, 
                               double lGP_error_variance, 
                               arma::mat Omega,
                               arma::mat invOmega,
                               arma::mat kk
                               )
{
  // p maximum length of the dependence in Omega
  
  int n = AA.n_cols ;
  int T = AA.n_rows ;
  int newT = kk.n_cols ; 
  
  arma::vec tmp_meanlatentGP(T, arma::fill::value(0.0)) ;


  // reconstruct the binary signal for each (i,t)
  // matrix T x n, binary_signal(t,i) = I( a(t,i) > 0 )
  arma::mat binary_signal(T, n, arma::fill::zeros) ;
  for(int i = 0; i < n; i++) {
    for(int t = 0; t < T; t++) {
      if(AA(t,i) > 0.0) { binary_signal(t,i) = 1.0 ; }
    }
  }

  /*
   * DATA AUGMENTATION z_it | ... ~ truncNorm
   */
  arma::mat z_tilde(T, n) ;
  for(int i = 0; i < n; i++) {
    for(int t = 0; t < T; t++) {
      if(binary_signal(t,i) > 0.5) {
        z_tilde(t, i) = r_truncnorm(latent_signal(t, cl_s(i)), 1, 0.0, 999.0) ; // truncated on (0, +infty)
      } else {
        z_tilde(t, i) = r_truncnorm(latent_signal(t, cl_s(i)), 1, -999.0, 0.0) ; // truncated on (-infty, 0)
      }
    }
  }
  
  
  /*
   * CLUSTER-SPECIFIC GP REALIZATIONS
   */
  for(int cli = 0; cli < max_cl; cli++)
  {
    arma::uvec idcl = arma::find(cl_s == cli) ;
    double n_k = idcl.n_elem ;

    if( n_k > 0 ) {
      arma::mat z_tilde_cl = z_tilde.cols(idcl) ;
      arma::vec mean_z_tilde_cl = arma::sum(z_tilde_cl, 1) / (n_k*1.0) ; // mean_{i} z_{it} vector length T
      
      arma::vec predicted_mean(T) ;
      arma::vec est_latent_signal(T) ;
      
      for(int j = 0 ; j < newT ; j++){
        arma::vec new_mean = kk.col(j).t() * invOmega * mean_z_tilde_cl ; 
        arma::vec new_var = (lGP_sigma2_omega + lGP_error_variance) - kk.col(j).t() * invOmega * kk.col(j) ; 
        
        est_latent_signal(j) = R::rnorm(new_mean(0), sqrt(new_var(0))) ; 
        predicted_mean(j) = new_mean(0) ; 
        
        latent_signal.col(cli) = predicted_mean ;
      }
      
      

    } else {
      arma::mat tmpN = rmvnorm(1, tmp_meanlatentGP, Omega) ;
      latent_signal.col( cli ) =  tmpN.row(0).t() ;
    }
  }
  
  return(latent_signal) ;
}


