#-------------------# #-------------------# #-------------------#
#                   PERFORM SPIKE DETECTION
#-------------------# #-------------------# #-------------------#

# These scripts replicate the comparison with the two-step approach reported in Section S4.3 of the Supplementary Material.
# Specifically, this script uses the L0 optimization method of Jewell and Witten (2020) to deconvolve the traces.

#  _________________________________
#  YOU CAN AVOID RUNNING THIS SCRIPT 
#  ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
# In the Google Drive folder (https://drive.google.com/drive/folders/1-1xf57mZBc1usA-iCZGkp4KPF8oX_5mV?usp=sharing)
# the RDS files produced by this script are available in the "Code/03_simulation_study/02_comparison_two-step_synthetic_data/results" folder.
# You can download these files and copy them in the corresponding folder of this repository.
# Notice that these files are not necessary for replicating the results of the paper.

replications_sim = 50

library(FastLZeroSpikeInference)

estimated_calcium = matrix(0, 100, 80)

for(n_sim in 1:replications_sim){
  nameopen = paste0("03_simulation_study/01_sensitivity_study/simulated_data/simulated_data", n_sim,".RDS")
  sim = readRDS(nameopen)
  
  TT = sim$TT
  N = sim$n

  dd = apply(t(sim$calcium), 2, diff)
  index_dd = dd > quantile(dd, .90)
  
  sigma2_start = var(sim$calcium[sim$calcium <  quantile(dd, .70)])

  #-----------------------# #-----------------------# 
  #      apply JW model for spike detection         #
  #-----------------------# #-----------------------# 
  lambdas = seq(1,100, length.out = 300)
  gam = 0.9
  sd_trace = sqrt(sigma2_start)
  
  # select the parameter lambda so that the number of detected spikes smaller than 1 estimated sd is zero
  fit = sapply(lambdas,
               function(lam) list(FastLZeroSpikeInference::estimate_spikes(dat = sim$calcium[1,],
                                                                           gam = gam, lambda = lam,
                                                                           estimate_calcium = T) ))

  mean_spike = (lapply(fit, function(x) 
    sum((x$estimated_calcium[x$spikes] - x$estimated_calcium[x$spikes-1])<sd_trace)  ))

  lambda = lambdas[which(unlist(mean_spike)==0)[1]]

  fit = list()
  for(i in 1:nrow(sim$calcium))
  {
    fitL0 = FastLZeroSpikeInference::estimate_spikes(dat = sim$calcium[i,],
                                                     gam = gam, lambda = lambda,
                                                     estimate_calcium = T)
    fit[[i]] = fitL0$spikes
    estimated_calcium[i,] = fitL0$estimated_calcium
  }
 
  estimated_signal = matrix(0, N, TT)
  for(i in 1:N){
    estimated_signal[i,fit[[i]]] = 1
  }
  
  
  
  name = paste0("03_simulation_study/02_comparison_two-step_synthetic_data/results/estimates_JW_simu", n_sim, ".RDS")
  saveRDS(estimated_signal, file = name)
  
  name = paste0("03_simulation_study/02_comparison_two-step_synthetic_data/results/estimated_calcium_JW_simu", n_sim, ".RDS")
  saveRDS(estimated_calcium, file = name)
  
  rm(estimated_signal)
  gc()

}
