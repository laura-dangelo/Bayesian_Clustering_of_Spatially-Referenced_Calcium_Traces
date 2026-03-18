#include "00_MH_gamma.h"



// prior on gamma
/*
 * gamma ~ Beta(hyp_gamma1, hyp_gamma2)
 */
double logprior_gamma(double gamma, double hyp_gamma1, double hyp_gamma2) 
{
  double out ;
  out = R::dbeta( gamma, hyp_gamma1, hyp_gamma2, true ) ;
  return(out) ;
}

// log-posterior (MH step on gamma)
double logpost_gamma(arma::mat Ca, arma::mat AA, 
                     double gamma, double tau2, 
                     double hyp_gamma1, double hyp_gamma2)
{
  int n = Ca.n_cols ;
  int T = Ca.n_rows -1 ;
  
  double out = 0 ;
  double llik = 0 ;
  for(int i  = 0; i < n; i++) {
    for(int t = 0; t < T; t++) { 
      llik = llik + R::dnorm(Ca(t+1,i), gamma*Ca(t,i) + AA(t,i), std::sqrt(tau2), true) ;
    }
  }
  out = llik + logprior_gamma(gamma, hyp_gamma1, hyp_gamma2) ;
  return(out) ;
}

