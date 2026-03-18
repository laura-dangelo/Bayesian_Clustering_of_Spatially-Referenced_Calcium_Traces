
replications_sim = 50


# devtools::install_github("jewellsean/FastLZeroSpikeInference")
library(FastLZeroSpikeInference)

estimated_calcium = matrix(0, 100, 80)

for(n_sim in 1:replications_sim){
  nameopen = paste0("10_simulations_sensitivity/simulated_data/simulated_data", n_sim,".RDS")
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
  # (matrix(unlist(mean_spike), dimnames = list(c(round(lambdas,4))) ))
  
  lambda = lambdas[which(unlist(mean_spike)==0)[1]]
  # lambda
  # fit = FastLZeroSpikeInference::estimate_spikes(dat = sim$calcium[1,],
  #                                                gam = gam,
  #                                                lambda = lambda,
  #                                                estimate_calcium = T)
  # 
  # fit$spikes # spike times
  # amplitudes = fit$estimated_calcium[fit$spikes] - fit$estimated_calcium[fit$spikes-1]
  # 
  # plot(sim$calcium[1,], type="l", lwd=2)
  # for(i in 1:length(fit$spikes)) segments(fit$spikes[i], 0, fit$spikes[i], amplitudes[i], col=2)


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
  
  
  
  name = paste0("11_simulations_comparison_two-step/results/estimates_JW_simu", n_sim, ".RDS")
  saveRDS(estimated_signal, file = name)
  
  name = paste0("11_simulations_comparison_two-step/results/estimated_calcium_JW_simu", n_sim, ".RDS")
  saveRDS(estimated_calcium, file = name)
  
  rm(estimated_signal)
  gc()

}
