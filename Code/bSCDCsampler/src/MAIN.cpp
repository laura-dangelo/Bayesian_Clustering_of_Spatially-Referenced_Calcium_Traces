
#include "00_cov_functions.h"
#include "00_sample_mvnorm.h"
#include "00_MH_gamma.h"

#include "01_cluster_neurons_PSB.h"
#include "01_cluster_amplitudes.h"
#include "01_sample_latent_signal.h"

// [[Rcpp::depends(RcppDist)]]
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>


// Gibbs sampler: main function
// [[Rcpp::export]]
Rcpp::List calcium_gibbsCPP_burnin(int Nrep, 
                                    int burnin,
                                    arma::mat Y, // matrix (T x n) y_ti
                                    arma::mat loc, // matrix (n x 2)
                                    double par_varcov_loc, // covariance between locations
                                    arma::vec cluster_signal_start, // starting point of the cluster on signal: vector length n
                                    int max_cl, // maximum number of  clusters 
                                    double b_start, double gamma_start, 
                                    double sigma2_start, // variance on observed calcium
                                    double tau2_start, // variance on latent calcium
                                    arma::mat cluster_amplitudes_start, // starting point of the cluster on amplitudes: matrix (n x T)
                                    arma::vec ak_start, // unique amplitudes: vector of length max_cl
                                    double varC0, // variance calcium at time 0
                                    double hyp_a1, double hyp_a2, // prior shape and rate on amplitudes
                                    double hyp_b1, double hyp_b2, // prior mean and variance on b (Gaussian)
                                    double hyp_sigma21, double hyp_sigma22, // prior shape and rate Gamma parameters on 1/sigma2 (precision)
                                    double hyp_tau21, double hyp_tau22, // prior shape and rate Gamma parameters on 1/tau2 (precision state equation)
                                    double hyp_gamma1, double hyp_gamma2, // prior shape1 and shape2 parameter Beta on gamma (AR decay param)
                                    double eps_gamma, // MH step size on gamma
                                    double step_AMH_gamma, // step adaptive MH on gamma
                                    double eps_a, // MH step size on amplitudes
                                    double step_AMH_a,
                                    double a_low, // threshold for minimum spike amplitude
                                    double lGP_tau2_omega, // length scale of the latent GP
                                    double lGP_sigma2_omega, //  variance of the latent GP
                                    double lGP_error_variance, // variance of the error of observations of the latent GP 
                                    double mu_alpha, double sigma2_alpha // parameters data augmentation probit SB
                                    ) 
{
  int n = Y.n_cols ;
  int T = Y.n_rows ;
  bool error = false ;
  
  /*
  * compute covariance matrices Sigma and Omega
  */
  // Sigma (varcov_loc): proximity matrix between neurons
  arma::mat varcov_loc(n,n) ;
  varcov_loc = compute_Sigma(loc, par_varcov_loc) ;
  
  // Omega: covariance matrix modeling the temporal dependence between spikes
  arma::mat mat_error_var(T, T, arma::fill::eye);
  mat_error_var =  mat_error_var * (lGP_error_variance/(n*1.0)) ;
  arma::mat Omega = arma::symmatu(compute_Omega(T, lGP_tau2_omega, lGP_sigma2_omega)) + mat_error_var ;
  arma::mat invOmega = inv_sympd(arma::symmatu(Omega)) ;
  
  arma::vec new_times = arma::linspace(0, T-1, T) ;
  arma::mat kk(T, new_times.n_elem) ;
  for(int j = 0 ; j < new_times.n_elem ; j++){
    kk.col(j) = compute_cor(new_times(j), T, lGP_tau2_omega, lGP_sigma2_omega) ;
  }
  
  
  // some tmp variables
  double tmpb ;

  // quantities for Kalman filter on calcium
  arma::mat out_Ca(T+1, n) ; // series of the calcium: I need the calcium at time 0, hence out_Ca[t,i] corresponds to y[t+1,i]
  out_Ca = arma::join_cols(arma::zeros(n).t(), Y) ; // (T+1) x n, first row is 0
  arma::vec filter_mean(T+1) ; // 0 1 ... n
  arma::vec filter_var(T+1) ; // 0 1 ... n
  arma::vec R(T+1) ; arma::vec a(T+1) ; // 0 1 ... n
  double back_mean; double back_var ;
  
  // probit stick-breaking
  arma::mat latent_alpha(n, max_cl) ;
  arma::mat inv_varcov_loc = inv_sympd(arma::symmatu(varcov_loc)) ;
  Rcpp::List out_pSB ;
  arma::vec tmp_mu(n) ; tmp_mu.zeros() ;
  latent_alpha = rmvnorm(max_cl, tmp_mu, varcov_loc).t() ;
  
  // sampling signal
  arma::mat logprobit_latentGP(T, max_cl, arma::fill::value( R::pnorm5(0.0, 0, 1, true, true) )) ;
  arma::mat out_statespace(T+1, max_cl, arma::fill::zeros) ;
  arma::vec log_SB_weights(max_cl-1) ; // stick-breaking probabilities of the cluster allocation of the non-zero amplitudes
  log_SB_weights = log_stick_breaking( Rcpp::rbeta(max_cl-1, 1.0, 1.0) ) ;
  
  
  /*
  * INITIALIZATION
  */
  // allocate matrices of MCMC output
  arma::mat out_b(n, Nrep - burnin) ; // mean of the series
  arma::vec out_gamma = arma::zeros(Nrep - burnin) ; // decay parameter
  arma::vec out_sigma2 = arma::zeros(Nrep - burnin) ; // variance on the output equation
  arma::vec out_tau2 = arma::zeros(Nrep - burnin) ; // variance on the state equation
  arma::mat out_cluster_signal(n, Nrep - burnin) ; // cluster allocation of the n neurons
  arma::cube out_cluster_amplitudes(T, n, Nrep - burnin, arma::fill::zeros) ; // cluster allocation of the amplitude parameters (0 if no spike, >0 if spike)
  arma::cube out_AA(T, n, Nrep - burnin, arma::fill::zeros) ; // matrix of the {0/a} signal indicator for each (i,t)
  arma::cube latent_signal(T, max_cl, Nrep - burnin, arma::fill::value(0.0)) ; // latent realizations of the GP underlying the binary signal
  arma::mat out_a(max_cl, Nrep - burnin, arma::fill::zeros) ; // unique values of the amplitudes 
  
  // initialize the tmp quantities of the burnin
  
  arma::vec tmp_b = arma::vec(n, arma::fill::value(b_start)) ;
  double tmp_gamma = gamma_start ;
  double tmp_sigma2 = sigma2_start ;
  double tmp_tau2 = tau2_start ;
  arma::vec tmp_cluster_signal = cluster_signal_start ;
  arma::vec tmp_a = ak_start ;
  arma::mat tmp_cluster_amplitudes = cluster_amplitudes_start ; 
  arma::mat tmp_AA(T, n, arma::fill::zeros) ;
  arma::mat tmp_latent_signal(T, max_cl, arma::fill::zeros) ;

  
  /*
  * Metropolis-Hastings steps
  */
  double ratio ;
  
  // MH on gamma (decay parameter)
  double gamma_old ; double gamma_new ;
  arma::vec acc_vec_MHgamma(Nrep + 1, arma::fill::zeros);
  
  // MH on unique amplitudes a*_k
  Rcpp::List out_DPmix_a ;
  arma::mat acc_mat_MHa(Nrep + 1, max_cl, arma::fill::zeros) ;
  arma::vec eps_a_v(max_cl); eps_a_v.fill(eps_a);
  
  
  //// progress bar ////
  bool display_progress = true ;
  Progress p2(Nrep, display_progress) ;
  
  /*
  * START MCMC
  */
  for(int iter = 0; iter < Nrep ; iter++)
  {
    if( Progress::check_abort() ) { return -1.0 ; }
    p2.increment();
    
    /*
    * Reconstruct the matrix s(i,t)*a(i,t) for all i,t
    */
    tmp_AA.zeros() ;
    for(int i = 0; i < n; i++) {
     for(int t = 0; t < T; t++) {
       if(tmp_cluster_amplitudes(t, i) > 0) {
         tmp_AA(t, i) = tmp_a( tmp_cluster_amplitudes(t, i) ) ;
       }
     }
    }
      
     
     
    /*
    * Sample calcium level c(i,t)
    * Kalman filter + backward sampling
    */
    for(int i = 0; i < n; i++)
    {
      a(0) = 0 ; // a0
      R(0) = varC0 ; // R0
      filter_mean(0) = 0 ;
      filter_var(0) = varC0 ;
      for(int t = 1; t < T + 1 ; t++) {
        a(t) = tmp_gamma * filter_mean(t-1) + tmp_AA(t-1, i);
        R(t) = pow(tmp_gamma,2) * filter_var(t-1) + tmp_tau2 ;
           
        filter_mean(t) = a(t) + R(t) / (R(t) + tmp_sigma2 ) * (Y(t-1, i) - tmp_b(i) - a(t)) ;
        filter_var(t) = tmp_sigma2 * R(t) / (R(t) + tmp_sigma2 ) ;
      }
      out_Ca(T, i) = R::rnorm(filter_mean(T), sqrt(filter_var(T))) ;
         
      for(int t = T-1; t > -1; t--)
      {
        back_mean = filter_mean(t) + tmp_gamma * filter_var(t) / R(t+1) * (out_Ca(t+1, i) - a(t+1)) ;
        back_var = filter_var(t) - pow(tmp_gamma * filter_var(t), 2) / R(t+1) ;
           
        out_Ca(t, i) = R::rnorm(back_mean, sqrt(back_var)) ;
      }
    }
      
    /*
    * Sample b = [b_1, ..., b_n]
    * semi-conjugate
    */
    arma::mat z(T,n) ; // z(i,t) = y(i,t) - Ca(i,t) difference between observed level and estimated calcium
    for(int i = 0; i < n; i++)
    {
      for(int t = 0; t < T; t++) { z(t,i) = Y(t,i) - out_Ca(t+1, i) ; }
      tmpb = Y.n_rows * hyp_b2 + tmp_sigma2 ;
      tmp_b(i) = R::rnorm( (hyp_b2 * arma::accu(z.col(i)) +  tmp_sigma2 * hyp_b1)/tmpb , // new mean
                                   std::sqrt( (hyp_b2 * tmp_sigma2) / tmpb ) // new s.d.
      ) ;
    }
      
      
    /*
    * Sample sigma2
    */
    arma::mat sq(T,n) ; 
    for(int t = 0; t < T; t++) { 
      // compute squares of (y_it - b_it - c_it )
      for(int i = 0; i < n; i++) {  sq(t,i) = (z(t,i) - tmp_b(i)) * (z(t,i) - tmp_b(i)) ; }
    }
    // parameterization of the gamma in rcpp with shape and rate (not scale)
    tmp_sigma2 = 1.0 / R::rgamma(hyp_sigma21 + (n*T*1.0)/2.0, 
                                 1.0 /(hyp_sigma22 + 0.5 * arma::accu(sq)) ) ; 
      
  
    /*
    * Sample tau2
    */
    arma::mat sq2(T,n) ;
    for(int t = 0; t < T; t++) {
      for(int i = 0; i < n; i++) {
        // compute squares of (c_it - gamma c_{i,t-1} - a_it)
        sq2(t,i) = (out_Ca(t+1,i) - tmp_gamma * out_Ca(t,i) - tmp_AA(t,i)) *
          (out_Ca(t+1,i) - tmp_gamma * out_Ca(t,i) - tmp_AA(t,i)) ;
      }
    }
    tmp_tau2 = 1.0 / R::rgamma(hyp_tau21 + (n*T*1.0)/2.0, 1.0/(hyp_tau22 + 0.5 * arma::accu(sq2)) ) ;
      
    
      
    /*
    * Sample gamma
    * MH step with random walk on logit transform
    */
    gamma_old = tmp_gamma ;
    double eta_old = log(gamma_old/(1.0 - gamma_old)) ;
    double eta_new = eta_old + R::rnorm(0.0, eps_gamma) ;
    gamma_new = exp(eta_new)/(1.0 + exp(eta_new)) ;
    ratio = logpost_gamma(out_Ca, 
                          tmp_AA, 
                          gamma_new, tmp_tau2,
                          hyp_gamma1, hyp_gamma2) -
            logpost_gamma(out_Ca, 
                          tmp_AA,
                          gamma_old, tmp_tau2,
                          hyp_gamma1, hyp_gamma2) +
            eta_new-eta_old - 2.0 * (log(1.0 + exp(eta_new))) + 2.0 * (log(1.0 + exp(eta_old))) ;
    if(R::runif(0, 1) < exp(ratio)) { 
      gamma_old = gamma_new ; 
      acc_vec_MHgamma(iter) = 1 ;
    }
    tmp_gamma = gamma_old ;
      
    // Adaptive MH
    double delta_n = std::min(step_AMH_gamma, 1.0/sqrt(iter)) ;
    if((iter % 50 == 0) & (iter > 49)){
      arma::vec batch = acc_vec_MHgamma(arma::span(iter-50, iter-1)) ;
      double prop_batch = mean(batch) ;
      if(prop_batch > .44) {
        eps_gamma =  exp(log(eps_gamma)+delta_n) ;
      } else {
        eps_gamma =  exp(log(eps_gamma)-delta_n) ;
      }
    }
      
      
      
    /*
    * Sample the neurons' cluster allocation using location-dependent probit SB process 
    */
    // Part 1: sample the cluster allocation variables
    out_pSB = cluster_neurons(tmp_AA,
                              logprobit_latentGP,
                              tmp_cluster_signal, max_cl,
                              inv_varcov_loc,
                              latent_alpha,
                              mu_alpha,
                              sigma2_alpha) ;
    
    if(out_pSB["error"]) { Rcpp::Rcout << "error pSB" ; error = true ; iter = Nrep; }
    
    arma::vec out_cluster_signalfun = out_pSB["cluster_signal"] ;
    tmp_cluster_signal = out_cluster_signalfun ;
    arma::mat out_lafun = out_pSB["latent_alpha"] ;
    latent_alpha = out_lafun ;
    
    
    // Part 2: update the atoms
    // atoms here are the latent continuous signal for each cluster
    // realizations of a GP indexed by time, with covariance function Omega
    out_statespace = sample_latent_signal(tmp_AA,
                                          tmp_latent_signal,
                                          tmp_cluster_signal,
                                          max_cl,
                                          lGP_tau2_omega, 
                                          lGP_sigma2_omega, 
                                          lGP_error_variance, 
                                          Omega,
                                          invOmega,
                                          kk) ;
    tmp_latent_signal = out_statespace ;
    
    // Part 3: compute the probit transformation of the latent GP signal
    for(int k = 0; k < max_cl; k++) {
      for(int t = 0; t < T; t++) {
        logprobit_latentGP(t,k) = R::pnorm( tmp_latent_signal(t, k) , 0, 1, true, true ) ;
      }
    }
    
  
    
    /*
     * Sample amplitudes
     */
    out_DPmix_a = cluster_amplitudes(out_Ca,
                                     tmp_gamma,
                                     tmp_latent_signal,
                                     tmp_cluster_signal,
                                     tmp_tau2,
                                     log_SB_weights,
                                     tmp_a,
                                     tmp_cluster_amplitudes,
                                     max_cl,
                                     1.0,
                                     hyp_a1, hyp_a2,
                                     eps_a_v, 
                                     step_AMH_a,
                                     iter,
                                     acc_mat_MHa,
                                     a_low
                                     ) ;
    
    if(out_DPmix_a["error"]) {  Rcpp::Rcout << "error amplitudes" ; error = true ; iter = Nrep; }
    arma::vec out_eps_a_v = out_DPmix_a["eps_a_vec"] ;
    eps_a_v = out_eps_a_v;
    arma::rowvec acc_it_a = out_DPmix_a["acc_a"] ;
    acc_mat_MHa.row(iter+1) = acc_it_a; 
    
    arma::mat output_cluster_amplitudes = out_DPmix_a["cluster_amplitudes"] ;
    tmp_cluster_amplitudes = output_cluster_amplitudes ;
    
    
    arma::vec output_a = out_DPmix_a["a"] ;
    tmp_a = output_a ;
    
    arma::vec output_prob = out_DPmix_a["log_SB_weights"] ;
    log_SB_weights = output_prob ;
    
  
  
    
    if(iter >= burnin) {
      out_AA.slice(iter - burnin) = tmp_AA ;
      out_b.col(iter - burnin) = tmp_b ;
      out_sigma2(iter - burnin) = tmp_sigma2 ;
      out_tau2(iter - burnin) = tmp_tau2 ;
      out_gamma(iter - burnin) = tmp_gamma ;
      latent_signal.slice(iter - burnin) = tmp_latent_signal ;
      out_cluster_amplitudes.slice(iter - burnin) = tmp_cluster_amplitudes ;
      out_a.col(iter - burnin) = tmp_a ;
      out_cluster_signal.col( iter - burnin ) = tmp_cluster_signal ;
    }
    
    
    // END GIBBS  
  }

  return Rcpp::List::create(Rcpp::Named("latent_signal") = latent_signal,
                            Rcpp::Named("AA") = out_AA,
                            Rcpp::Named("amplitudes") = out_a,
                            Rcpp::Named("cluster_signal") = out_cluster_signal,
                            Rcpp::Named("cluster_amplitudes") = out_cluster_amplitudes,
                            Rcpp::Named("Ca") = out_Ca,
                            Rcpp::Named("b") = out_b,
                            Rcpp::Named("gamma") = out_gamma,
                            Rcpp::Named("sigma2") = out_sigma2,
                            Rcpp::Named("tau2") = out_tau2,
                            Rcpp::Named("error") = error);
}
