#include "00_cov_functions.h"

// covariance matrix on the Gaussian process
// compute Omega given the length of the series TT
// the length-scale lGP_tau2_omega, and the variance of the fluctuations lGP_sigma2_omega
// [[Rcpp::export]]
arma::mat compute_Omega(int TT, double lGP_tau2_omega, double lGP_sigma2_omega)
{
  arma::mat out(TT,TT) ;
  arma::vec times = arma::linspace(0, TT-1, TT);
  arma::vec tmp(TT) ;
  
  for(int t = 0; t < TT; t++) { tmp(t) = lGP_sigma2_omega * exp(- ( times(t)*times(t) ) / (2.0 * lGP_tau2_omega)  ) ; }
  for(int t = 0; t < TT; t++)
  {
    for(int j = 0; j < TT; j++)
    {
      out(t,j) = tmp( abs(t-j) ) ;
    }
  }
  return(out) ;
}

// compute Sigma(theta) given theta and the locations
// proximity matrix 

// [[Rcpp::export]]
arma::mat compute_Sigma(arma::mat loc, double theta)
{
  int n = loc.n_rows ;
  arma::mat out(n, n, arma::fill::ones) ;
  for(int i = 1; i < n; i++){
    for(int j = 0; j < i; j++){
      out(i,j) = out(j,i) = exp(- (.5 / theta) * ( (loc(i,0)-loc(j,0))*(loc(i,0)-loc(j,0)) + 
        (loc(i,1)-loc(j,1)) * (loc(i,1)-loc(j,1)) ) ) ;
    }
  }
  return(out) ;
}


arma::vec compute_cor(double t_new, int T, double lGP_tau2_omega, double lGP_sigma2_omega)
{
  arma::vec out(T) ;
  arma::vec times = arma::linspace(0, T-1, T);
  
  for(int j = 0; j < T; j++) {
    double diff = t_new - times(j)*1.0  ;
    out(j) = lGP_sigma2_omega * exp(- (diff * diff) / (2.0 * lGP_tau2_omega)  ) ;
  }
  
  return(out) ;
}