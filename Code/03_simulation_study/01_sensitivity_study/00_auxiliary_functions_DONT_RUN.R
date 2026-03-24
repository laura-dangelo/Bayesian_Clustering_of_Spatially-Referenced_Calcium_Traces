
#---------#    AUXILIARY FUNCTIONS     #---------#
#---------#   don't run this script    #---------#


# squared exponential kernel covariance function with unit variance
rbf_kernel = function(theta, x1, x2)
{
  arg = - (.5 / theta) * sum((x1-x2)^2)
  return(exp(arg))
}

# simulate a realization from a GP of parameters
# tau2: lengthscale
# sigma2: overall variance (amplitude)
# mu: mean (constant)
varcov_GP = function(len, par_tau2, par_sigma2)
{
  x_seq = 0:(len-1)
  x_exp = as.matrix(expand.grid(x_seq, x_seq))
  par_sigma2 * matrix(apply(x_exp, 1, function(x) rbf_kernel(theta = par_tau2, x1 = x[1], x2 = x[2])), len, len)
}

sim_GP = function(mean, varcov)
{
  len = nrow(varcov)
  out = rmvnorm(1, mean, varcov)
  return(out)
}



# generate the sequences of binary signal
gen_signal = function(n, TT, cluster,
                      meanGP,
                      par_tau2, 
                      par_sigma2) 
{
  trunc = length(unique(cluster))
  atoms_GP = matrix(NA, trunc, TT) 
  Omega = varcov_GP(TT, par_tau2 = par_tau2, par_sigma2)
  
  ## generate atoms GP
  for(k in 1:trunc) atoms_GP[k,] = sim_GP(meanGP[k,], Omega)
  
  ## generate signal
  signalk = matrix(NA, n, TT)
  for(i in 1:n) {
    signalk[i,] = rbinom(TT, size = 1, prob = pnorm(atoms_GP[cluster[i],]) ) 
  }
  # signalk[,c(1:3)] = 0
  
  return(list(atoms_GP = atoms_GP, signal = signalk))
}


## signal into trace
signal_to_trace = function(signal, gamma, tau2, sigma2)
{
  out = numeric(length(signal)+1)
  out[1] = rnorm(1, 0, sqrt(tau2))
  for(i in 2:(length(out)))
  {
    out[i] = gamma * out[i-1] + signal[i-1] + rnorm(1, 0, sqrt(tau2))
  }
  return(out[-1] + rnorm(length(signal), 0, sqrt(sigma2)))
}

# generate amplitudes of the spikes, generate the AR(1) process, add noise
gen_amplitudes_trace = function(signal,
                                ampl, # this is a vector
                                gamma, tau2, sigma2) 
{
  n = nrow(signal)
  TT = ncol(signal)
  clusters_amplitudes = sample(1:length(ampl), n*TT, replace = T)
  amplitudes = matrix( ampl[clusters_amplitudes], n, TT)
  
  sig_ampl = signal * amplitudes
  clus_ampl = signal * clusters_amplitudes
  
  ## generate calcium trace
  calcium = t(apply(sig_ampl, 1, function(x) signal_to_trace(x, gamma, tau2, sigma2))) 
  
  return(list(amplitudes = sig_ampl, clus_ampl = clus_ampl, trace = calcium))
}