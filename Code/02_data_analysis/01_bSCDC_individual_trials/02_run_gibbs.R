
#---------# #-----------# #---------# #---------# #-----------# #---------# 
#---------#                  RUN GIBBS SAMPLER                  #---------#
#---------# #-----------# #---------# #---------# #-----------# #---------# 

# This script runs the Gibbs sampler algorithm on all the time windows. This can be time consuming. 
# To allow you to try it on a smaller subset of data, the flag variable run_on_subset is now set to TRUE. 
# Under this configuration, the sampler is ran only on the time windows reported in the article (main and supplementary).
# If you want to further speed up computation, you can consider reducing the number of MCMC iteration, which is now
# set to quite a large value to ensure good convergence.

# If you wish instead to proceed with the full computation, set run_on_subset = FALSE.
# Notice that to produce the plots for the overall activity you will need the output RDS files for all the windows.
# The pre-computed RDS are available as Supplementary Material. 

install.packages("bSCDCsampler_0.0.1.tar.gz")
library(bSCDCsampler)

idx = readRDS("../Data/Time_windows/indices.RDS")
loc_neurons = readRDS("../Data/M3424F_loc_neurons.RDS")

# These variables define the number of MCMC iterations and burn-in
nrep = 30000
burnin = 25000

# set this variable to FALSE to run all the windows
run_on_subset = TRUE

if(run_on_subset) {
  idx_to_run = c(4,17,36,77)
} else {
  idx_to_run = idx
}


#------# #------# #------# #------# #------#
#------#   fixed hyperparameters    #------#
#------# #------# #------# #------# #------#

# Prior on the amplitudes
par1 = 10 # shape hyperparameter of the Gamma prior on a
par2 = 1 # rate hyperparameter of the Gamma prior on a
a_low = 1 # threshold 

# Prior on the PSBP (location-informed prior for clustering)
par_varcov_loc = 500 # theta: parameter of the PSBP governing the closeness between neurons
par_tau2 = 10 # lengthscale parameter of the latent GP
par_sigma2 = 4 # variance parameter of the latent GP

#------# #------# #------# #------# #------#



#------# #------# #------# #------# #------#
#------#      run Gibbs sampler     #------#
#------# #------# #------# #------# #------#

n_w = 0 
n_window = idx[n_w+1]
for(n_window in idx_to_run) {
  n_w = n_w+1
  cat(paste("Running window",n_w,"out of",length(idx_to_run)," - ID:",n_window,"\n"))
  
  filename = paste0("../Data/Time_windows/calcium_window", n_window, ".RDS")
  calcium = readRDS(filename)
  
  TT = nrow(calcium)
  N = ncol(calcium)
  trunc = 50 # maximum number of clusters
  
  set.seed(123*n_w) 
  
  
  # To ease convergence, we initialize the variances to reasonable values
  # Initialization of the variances can be delicate. We suggest the following, which provides stable and good results:
  dd = apply(calcium, 2, diff)
  index_dd = dd > quantile(dd, .95)
  sigma2_start = var(t(calcium)[t(calcium) <  quantile(dd, .8)])
  tau2_start = sigma2_start/100
  
  # We also provide reasonable starting points for the amplitudes of the spikes and the cluster assignments
  a_clus = kmeans( dd[index_dd], trunc-1)
  a_star = c(0, a_clus$centers)
  cl_a_start = matrix(0, TT-1, N)
  cl_a_start[index_dd] = a_clus$cluster 
  cl_a_start = rbind(rep(0,N), cl_a_start)
  cl_neurons = kmeans(t(calcium), 10)
  cl_s_start = cl_neurons$cluster -1
  
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
  filenameout = paste0("02_data_analysis/01_bSCDC_individual_trials/results/run_gibbs_window", n_window,"_alow_1.RDS")
  
  saveRDS(out, file = filenameout)
  
  cat(paste("Now removing old run...\n"))
  rm(out)
  rm(calcium)
  rm(a_clus)
  gc()
}

