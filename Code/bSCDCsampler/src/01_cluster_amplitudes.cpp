#include "01_cluster_amplitudes.h"


Rcpp::List cluster_amplitudes(arma::mat calcium,
                              double gamma, 
                              arma::mat latent_signal, 
                              arma::vec cl_s,
                              double tau2,
                              arma::vec log_SB_weights,
                              arma::vec unique_ampl,
                              arma::mat cl_a,
                              int max_cl,
                              double alpha,
                              double hyp_a1, double hyp_a2,
                              arma::vec eps_a_v, 
                              double step_AMH_a,
                              int iter,
                              arma::mat acc_mat_MHa,
                              double a_low
                              ) 
{
  int n = calcium.n_cols ;
  int T = calcium.n_rows -1 ;
  
  bool error = false ;
  
  // matrix of the latent calcium without the autoregressive behavior
  // cell (t,i) is Ca(t,i) - gamma * Ca(t-1,i)
  // calcium is initialized at time 0
  arma::mat calcium_noAR(T, n) ;
  for(int i = 0; i < n; i++) {
    calcium_noAR.col(i) = calcium(arma::span(1,T), i) - arma::vec(T, arma::fill::value(gamma)) % calcium(arma::span(0,T-1), i);
  }
  
  double logprob_spike ;
  arma::vec log_prob(max_cl) ; // log-posterior probability of cluster allocation
 
  arma::vec v_k(max_cl-1) ; // auxiliary beta r.v.'s
  arma::vec prob(max_cl, arma::fill::zeros) ;


  double old_a ; double new_a ; double ratio ;
  int a_k = 0 ; double b_k = 0.0 ; 
 
  /*
   * This is step (a) in the algorithm
   * sample cluster allocation
   */
  for(int i = 0; i < n; i++) {
    for(int t = 0; t < T; t++) {
      
      // posterior probability of not observing a spike 
      log_prob(0) = R::pnorm( -latent_signal(t, cl_s(i)), 0.0, 1.0, true, true ) + R::dnorm(calcium_noAR(t,i), 0.0, std::sqrt(tau2), true) ; 
      
      // logprob_spike is log( Phi(tilde{s}*_{k,t}) )
      // that is, the prior probability of observing a spike 
      logprob_spike = R::pnorm( latent_signal(t, cl_s(i)), 0.0, 1.0, true, true ) ; 

      // posterior probability of observing a spike of amplitude a*_k
      for(int k = 1; k < max_cl; k++) {
        log_prob(k) = logprob_spike + 
                      log_SB_weights(k-1) + 
                      R::dnorm(calcium_noAR(t,i), unique_ampl(k), std::sqrt(tau2), true) ;
      }
      prob = exp(log_prob - max(log_prob)) ;
    
      // sample new cluster allocation
      arma::vec cluster_id = arma::linspace(0, max_cl-1, max_cl) ;
      cl_a(t, i) = sample_cl(cluster_id, prob) ;
    }
  }
  
  
  arma::rowvec acc_vec_MHa(max_cl, arma::fill::zeros); 
  /*
   * This is step (b) in the algorithm
   * update cluster-specific parameter values a*_k
   */
  for(int k = 1; k < max_cl; k++)  // starts from 1 since k = 0 corresponds to the absence of a spike
  {
    
    arma::uvec ind_k = find( cl_a == k ) ;
    int nk = ind_k.n_elem ;
    if( nk > 0 )
    {

      // MH step: random walk with Gaussian transition density (tuning parameter is the s.d. eps_a_v)
      old_a = unique_ampl(k) ;
      double old_z = log(old_a-a_low);
      double new_z = R::rnorm(old_z, eps_a_v(k-1)) ;
      new_a = exp(new_z)+a_low;

      double llik_old = 0.0 ; double llik_new = 0.0 ;
      for(int j = 0; j < nk; j++) {
        llik_old = llik_old + R::dnorm( calcium_noAR(ind_k(j)), old_a, sqrt(tau2), true ) ;
        llik_new = llik_new + R::dnorm( calcium_noAR(ind_k(j)), new_a, sqrt(tau2), true ) ;
      }
      
      /* The ratio on the transformed parameter
       * ratio = posterior(new) / posterior(old)
       *       = (prior(new)*lik(new)) / (prior(old)*lik(old))
       *       = exp( logprior(new) + loglik(new) - logprior(old) - loglik(old) )
       */
      ratio =  R::dgamma(exp(new_z), hyp_a1, 1.0/hyp_a2, true) +
                   llik_new - llik_old -
               R::dgamma(exp(old_z), hyp_a1, 1.0/hyp_a2, true) +
                   new_z - old_z ;

      if(log(R::runif(0, 1)) < ratio) {
        old_a = new_a ; 
        acc_vec_MHa(k-1) = 1.0 ;
      }
      unique_ampl(k) = old_a ;

    } else {
      // this is the "else" part of condition "if( nk > 0 )"
      // sample from the prior
      unique_ampl(k) = R::rgamma(hyp_a1, 1.0/hyp_a2)+a_low ;
      acc_vec_MHa(k-1) = 1.0 ;
    }

    double delta_n = std::min(step_AMH_a, 1.0/sqrt(iter)) ;
    if((iter % 50 == 0) & (iter > 49)) {
      arma::vec batch = acc_mat_MHa(arma::span(iter-50,iter-1), k-1) ;
      double prop_batch = mean(batch) ;
      if(prop_batch > .44) {
        eps_a_v(k-1) =  exp(log(eps_a_v(k-1))+delta_n) ;
      } else {
        eps_a_v(k-1) =  exp(log(eps_a_v(k-1))-delta_n) ;
      }
    }
  }

  
  /*
   * This is step (c) in the algorithm
   * update the stick-breaking weights
   */
  for(int k = 1; k < max_cl; k++)
  {
    arma::uvec ind_k = find(cl_a == k) ;
    arma::uvec ind_mk = find(cl_a > k) ;
    
    a_k = 1 + ind_k.n_elem ;
    b_k = alpha + ind_mk.n_elem ;
    v_k(k-1) = R::rbeta(a_k, b_k) ;
  }
  log_SB_weights = log_stick_breaking( v_k ) ;
  
  
  return Rcpp::List::create(Rcpp::Named("a") = unique_ampl,
                            Rcpp::Named("cluster_amplitudes") = cl_a,
                            Rcpp::Named("log_SB_weights") = log_SB_weights,
                            Rcpp::Named("error") = error,
                            Rcpp::Named("eps_a_vec") = eps_a_v,
                            Rcpp::Named("acc_a") = acc_vec_MHa);
  
}
