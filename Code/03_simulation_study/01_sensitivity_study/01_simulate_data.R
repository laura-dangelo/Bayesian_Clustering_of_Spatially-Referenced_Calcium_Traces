
#-------------------# #-------------------# #-------------------#
##                   GENERATE SIMULATED DATA
#-------------------# #-------------------# #-------------------#

# These scripts replicate the sensitivity study reported in Section S4.2 of the Supplementary Material.
# Specifically, this script generates the synthetic data using the data generating mechanism outlined in Section S4.1.

#  _________________________________
#  YOU CAN AVOID RUNNING THIS SCRIPT 
#  ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
# the RDS files produced by this script are available in the "Code/03_simulation_study/01_sensitivity_study/simulated_data" folder.


library(mvtnorm)
source("03_simulation_study/01_sensitivity_study/00_auxiliary_functions_DONT_RUN.R")


# replications_sim: number of simulation replicates
replications_sim = 50

# n: number of neurons
# TT: length of the time series
n = 100
TT = 80 # average length of the windows on real data is 78



for(n_sim in 1:replications_sim)
{

  set.seed(n_sim)
  
  locations = rbind(
    rmvnorm(round(n/4), c(0, 200), diag(2)*1000),
    rmvnorm(n - round(n/3)- round(n/4), c(110, 220), diag(2)*500),
    rmvnorm(round(n/3), c(80, 150), diag(2)*700)
  )
  cluster_neurons = c( rep(1, round(n/4)), rep(2, n - round(n/3)- round(n/4)), rep(3, round(n/3)) )
  index_mixed = sample(1:n, n*(25/100), replace=F)
  for(i in 1:length(index_mixed)){
    cluster_neurons[index_mixed[i]] = sample( setdiff(1:3, cluster_neurons[index_mixed[i]]), 1 )
  }
  
  
  # parameters of the latent GP
  # length scale tau2: how wiggly are observations (amplitude w.r.t. t)
  # variance sigma2: amplitude of observations (amplitude w.r.t y)
  step_fun = matrix(-2, 3, TT)
  step_fun[1, c(1:4)] = -5
  step_fun[1, c(20:25)] = 3
  step_fun[1, c(49:53)] = 4
  
  step_fun[2, 45] = 3
  step_fun[2, 46:47] = 9
  step_fun[2, 48] = 1
  
  step_fun[3, 62:63] = 7
  step_fun[3, 65:66] = 2
  
  
  meanGP = t(apply(step_fun, 1, function(x) lm(x~poly(1:TT,7))$fitted.values ))*4

  par_tau2 = 5
  par_sigma2 = 1
  
  # generate signal
  sig1 = gen_signal(n, TT, cluster_neurons, meanGP = meanGP, par_tau2=par_tau2, par_sigma2 = par_sigma2)
  signal = sig1$signal
  atoms = pnorm(sig1$atoms_GP)
  
  # spike smplitudes
  unique_amplitudes = c(5,
                        6,6,6,
                        10,10,10,10,
                        14,
                        25) 
  
  # given the signal, generate the calcium traces 
  cal = gen_amplitudes_trace(signal, ampl = unique_amplitudes, gamma = 0.9, tau2 = 1, sigma2 = 1.5) 
  
  sim = list(n = n, TT = TT,
             calcium = cal$trace,
             loc = locations,
             cluster_neurons = cluster_neurons,
             signal = signal,
             amplitudes = cal$amplitudes,
             cluster_a = cal$clus_ampl,
             unique_amplitudes = unique_amplitudes,
             latent_prob = atoms,
             atoms_GP = sig1$atoms_GP
  )
  name = paste0("03_simulation_study/01_sensitivity_study/simulated_data/simulated_data", n_sim, ".RDS")

  saveRDS(sim, file = name)

}


