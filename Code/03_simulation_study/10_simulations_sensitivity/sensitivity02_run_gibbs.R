#-------------------# #-------------------# #-------------------#
#                   RUN GIBBS AND SAVE RESULTS
#-------------------# #-------------------# #-------------------#


library(salso)
library(corrplot)
library(bSCDCsampler)


#------# parameters of the Gamma prior on the amplitudes we are studying #------#

seq_alow = c(0, 0.5, 1) # threshold a_low
id_pars_gamma = c(1,2,3)
pars_gamma = matrix(c(3, 0.1, 4, 1, 10, 1), 3, 2, byrow = T) # parameters of the Gamma

parameter_comb = expand.grid(x = seq_alow, y = id_pars_gamma)
parameter_comb[,3] = sapply(parameter_comb[,2], function(x) pars_gamma[x,2])
parameter_comb[,2] = sapply(parameter_comb[,2], function(x) pars_gamma[x,1])
parameter_comb

#------# #------# #------# #------# #------# #------# #------# #------# #------# 

par_tau2 = 1
par_sigma2 = 2.5

replications_sim = 50

for(gg in 1:nrow(parameter_comb)){
  for(n_sim in 1:replications_sim){
    nameopen = paste0("10_simulations_sensitivity/simulated_data/simulated_data", n_sim,".RDS")
    sim = readRDS(nameopen)
    

    TT = sim$TT
    N = sim$n
    trunc = 30 # nrow(sim$atoms_GP)
    
    set.seed(123*n_sim) 
    
    dd = apply(t(sim$calcium), 2, diff)
    index_dd = dd > quantile(dd, .95)
    
    sigma2_start = var(sim$calcium[sim$calcium <  quantile(dd, .80)])
    tau2_start = sigma2_start/100
    
    a_clus = kmeans( dd[index_dd], trunc-1)
    a_star = c(0, a_clus$centers)
    
    cl_a_start = matrix(0, TT-1, N)
    cl_a_start[index_dd] = a_clus$cluster  
    cl_a_start = rbind(rep(0,N), cl_a_start)
    
    cl_neurons = kmeans((sim$calcium), 10)
    cl_s_start = cl_neurons$cluster -1
    
  
    nrep = 15000
    burnin = floor(nrep/5*2)
    
    par_varcov_loc = 500
    # corrplot(compute_Sigma_R(sim$loc[1:N,], par_varcov_loc), method = "color")
  
    
    tryCatch(
      {
        t1 <- Sys.time()
        out = calcium_gibbs_burnin(Nrep = nrep, 
                                   burnin = burnin,
                                   Y = t(sim$calcium)[1:TT, 1:N],  # matrix (T x n)
                                   loc = sim$loc[1:N,], # neurons' locations in the hippocampus
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
                                   hyp_a1 = parameter_comb[gg,2], hyp_a2 = parameter_comb[gg,3], # hyperparameters gamma prior on amplitudes
                                   hyp_b1 = 0, hyp_b2 = 1, # hyperparameters of the normal on the mean
                                   hyp_sigma21 = 3, hyp_sigma22 = 2, # hyperparameters of the invGamma on calcium variance
                                   hyp_tau21 = 10, hyp_tau22 = 1, # hyperparameters of the invGamma on state equation
                                   hyp_gamma1 = 3, hyp_gamma2 = 2, # hyperparameters of the beta on decay
                                   eps_gamma = 0.003, step_AMH_gamma = 0.005, # step MCMC on decay
                                   eps_a = 0.01, step_AMH_a = 0.01, # step MCMC on amplitude
                                   a_low = parameter_comb[gg,1], # shift of gamma on amplitudes
                                   lGP_tau2_omega = par_tau2, # lengthscale parameter of the latent GP
                                   lGP_sigma2_omega = par_sigma2, # variance parameter of the latent GP
                                   lGP_error_variance = 1, # output variance realizations of the latent GP
                                   mu_alpha = -3, sigma2_alpha = .5 # hyperparameters data augmentation probit SB
                                   )
        # str(out)
        t2 <- Sys.time()
        
        name = paste0("10_simulations_sensitivity/results/run_gibbs_alow_", 
                      parameter_comb[gg,1], "_", parameter_comb[gg,2], "_", parameter_comb[gg,3], 
                      "_simu", n_sim, ".RDS")
        saveRDS(out, file = name)
        name2 = paste0("10_simulations_sensitivity/results/time_gibbs_alow_", 
                      parameter_comb[gg,1], "_", parameter_comb[gg,2], "_", parameter_comb[gg,3], 
                      "_simu", n_sim, ".RDS")
        saveRDS(t2-t1, file = name2)
        rm(out)
        gc()
        
      }, error = function(e) {
        nameerror = paste0("results/ERROR_par", parameter_comb[gg,], "_run_", n_sim, ".Rdata")
        errorfile = c(0)
        save(errorfile, file = nameerror)
        print(nameerror)
      } )
  }
}

