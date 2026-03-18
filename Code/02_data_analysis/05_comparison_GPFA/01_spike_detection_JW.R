#---------# #-----------# #---------# #---------# #-----------# #---------# 
#---------#                DATA PREPARATION GPFA                #---------#
#---------# #-----------# #---------# #---------# #-----------# #---------# 
# deconvolution via L0 optimization.
# GPFA assumes a Negative Binomial likelihood. To return the number of spikes instead of just the presence/absence,
# we estimate the unitary spike and count how many spikes there are in a calcium increase.

library(FastLZeroSpikeInference)

calcium = readRDS("../Data/M3424F_calcium.RDS")

dd = diff(t(calcium))
sigma2_start = var(calcium[calcium <  quantile(dd, .95)])


lambdas = seq(30,350, length.out = 300)
gam = 0.95
sd_trace = sqrt(sigma2_start)

chosen_lambdas = numeric(325)
rm(dd)



for(i in 1:325){
  calcium = readRDS("../Data/M3424F_calcium.RDS")
  
  calcium = calcium[i,]

  fit = sapply(lambdas,
               function(lam) list(FastLZeroSpikeInference::estimate_spikes(dat = calcium,
                                                                           gam = gam, lambda = lam,
                                                                           estimate_calcium = T) ) )
  mean_spike = (lapply(fit, function(x) 
    sum((x$estimated_calcium[x$spikes] - x$estimated_calcium[x$spikes-1])<sd_trace)  ))
  summary(unlist(mean_spike))
  
  if(min(unlist(mean_spike))==00){
    lambda = lambdas[which(unlist(mean_spike)==0)[1]]
  } else {
    lambda = max(lambdas)
  }
  chosen_lambdas[i] = lambda
  
  fitL0 = FastLZeroSpikeInference::estimate_spikes(dat = calcium,
                                                   gam = gam, lambda = lambda,
                                                   estimate_calcium = T)
  
  estimated_amplitudes = numeric(length(fitL0$spikes))
  for(k in 1:length(estimated_amplitudes)){
    estimated_amplitudes[k] = fitL0$estimated_calcium[fitL0$spikes[k]] - fitL0$estimated_calcium[fitL0$spikes[k]-1]
  }
  number_spikes = round(estimated_amplitudes/sd_trace)
  
  
  name = paste0("02_data_analysis/05_comparison_GPFA/results/list_spikes_JW_neuron", i, ".RDS")
  saveRDS(fitL0$spikes, file = name)
  
  name = paste0("02_data_analysis/05_comparison_GPFA/results/estimated_amplitudes_JW_neuron", i, ".RDS")
  saveRDS(estimated_amplitudes, file = name)
  
  name = paste0("02_data_analysis/05_comparison_GPFA/results/number_spikes_JW_neuron", i, ".RDS")
  saveRDS(number_spikes, file = name)
  
  print(i)
  
  rm(calcium)
  rm(fit)
  rm(fitL0)
  rm(mean_spike)
}

saveRDS(chosen_lambdas, file = "02_data_analysis/05_comparison_GPFA/results/chosen_lambdas.RDS")

