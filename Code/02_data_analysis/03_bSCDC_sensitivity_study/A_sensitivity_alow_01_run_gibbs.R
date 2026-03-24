
#---------# #-----------# #---------# #---------# #-----------# #---------# #---------# 
#---------#         RUN GIBBS SAMPLER FOR SENSITIVITY ANALYSIS A_LOW        #---------#
#---------# #-----------# #---------# #---------# #-----------# #---------# #---------# 

# First part of the sensitivity analysis reported in Section S2.1 of the Supplementary Material. 
# This script runs the Gibbs sampler algorithm on window 17 for varying values of the threshold a_low.

# We study three different values for the threshold: 0, 0.5, 1
# Here, we are omitting the computation for the value 1 since it has already been computed in the data_analysis scripts

#  _________________________________
#  YOU CAN AVOID RUNNING THIS SCRIPT 
#  ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
# In the Google Drive folder (https://drive.google.com/drive/folders/1-1xf57mZBc1usA-iCZGkp4KPF8oX_5mV?usp=sharing)
# the RDS files produced by this script are available in the "Code/02_data_analysis/03_bSCDC_sensitivity_study/results" folder.
# You can download these files and copy them in the corresponding folder of this repository.


library(bSCDCsampler)

idx = readRDS("../Data/Time_windows/indices.RDS")
loc_neurons = readRDS("../Data/M3424F_loc_neurons.RDS")


par1 = 10
par2 = 1

a_low_seq = c(0, 0.5)

par_tau2 = 10
par_sigma2 = 4

n_window = 17  ### select window

for(a_low in a_low_seq) {

  filename = paste0("../Data/Time_windows/calcium_window", n_window, ".RDS")
  calcium = readRDS(filename)
  
  TT = nrow(calcium)
  N = ncol(calcium)
  trunc = 50 
  
  set.seed(123*(a_low+1)) 
  
  dd = apply(calcium, 2, diff)
  index_dd = dd > quantile(dd, .95)
  
  sigma2_start = var(t(calcium)[t(calcium) <  quantile(dd, .8)])
  tau2_start = sigma2_start/100
  
  a_clus = kmeans( dd[index_dd], trunc-1)
  a_star = c(0, a_clus$centers)
  
  cl_a_start = matrix(0, TT-1, N)
  cl_a_start[index_dd] = a_clus$cluster 
  cl_a_start = rbind(rep(0,N), cl_a_start)
  
  cl_neurons = kmeans(t(calcium), 10)
  cl_s_start = cl_neurons$cluster -1
  
  par_varcov_loc = 500
  nrep = 30000
  burnin = 25000
  
  out = calcium_gibbs_burnin(Nrep = nrep, 
                             burnin = burnin,
                             Y = as.matrix(calcium),  # matrix (T x n)
                             loc = loc_neurons, # neurons' locations in the hippocampus
                             par_varcov_loc = par_varcov_loc, # parameter of the probit SB governing the closeness between neurons
                             cluster_signal_start = cl_s_start, # initialization of the cluster of neurons
                             max_cl = trunc, # maximum number of mixture components (truncation parameter)
                             b_start = 0,  # initialization of the calcium mean
                             gamma_start = 0.9, # initialization of calcium decay
                             sigma2_start = sigma2_start, # initialization of the calcium variance
                             tau2_start = tau2_start, # initialization variance of the state equation calcium
                             cluster_amplitudes_start = cl_a_start, # initialization of the cluster of amplitudes
                             ak_start = a_star, # initialization of the unique values of the amplitudes
                             varC0 = 0.01, # variance of calcium at t = 0 (for Kalman filter)
                             hyp_a1 = par1, hyp_a2 = par2, # hyperparameters gamma prior on amplitudes
                             hyp_b1 = 0, hyp_b2 = 0.1, # hyperparameters of the normal on the mean
                             hyp_sigma21 = 3, hyp_sigma22 = 2, # hyperparameters of the invGamma on calcium variance
                             hyp_tau21 = 10, hyp_tau22 = 1, # hyperparameters of the invGamma on state equation
                             hyp_gamma1 = 3, hyp_gamma2 = 2, # hyperparameters of the beta on decay
                             eps_gamma = 0.005, step_AMH_gamma = 0.005, # step MCMC on decay
                             eps_a = 0.01, step_AMH_a = 0.01, # step MCMC on amplitude
                             a_low = a_low, # shift of gamma on amplitudes
                             lGP_tau2_omega = par_tau2, # lengthscale parameter of the latent GP
                             lGP_sigma2_omega = par_sigma2, # variance parameter of the latent GP
                             lGP_error_variance = 1, # output variance realizations of the latent GP
                             mu_alpha = -3, sigma2_alpha = .5 # hyperparameters data augmentation probit SB
  )
  
  cat(paste("Now saving...\n"))
  filenameout = paste0("02_data_analysis/03_bSCDC_sensitivity_study/results/run_gibbs_window", n_window,"_alow_", a_low,".RDS")
  
  saveRDS(out, file = filenameout)
  
  cat(paste("Now removing old run...\n"))
  rm(out)
  rm(calcium)
  rm(a_clus)
  gc()
}

