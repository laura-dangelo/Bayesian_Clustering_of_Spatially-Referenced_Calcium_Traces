
#---------# #-----------# #---------# #---------# #-----------# #---------# #---------# 
#---------#              STEP 1: DECONVOLUTION OF CALCIUM TRACES            #---------#
#---------# #-----------# #---------# #---------# #-----------# #---------# #---------# 
#
# This script performs deconvolution of the fluorescent calcium traces using the L0 optimization 
# algorithm of Jewell and Witten (2020), contained in the R package FastLZeroSpikeInference.
#
#  _________________________________
#  YOU CAN AVOID RUNNING THIS SCRIPT 
#  ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
# In the Google Drive folder (https://drive.google.com/drive/folders/1-1xf57mZBc1usA-iCZGkp4KPF8oX_5mV?usp=sharing)
# the RDS files produced by this script are available in the "Code/02_data_analysis/04_comparison_two-step/results" folder.
# You can download these files and copy them in the corresponding folder of this repository.



# devtools::install_github("jewellsean/FastLZeroSpikeInference")
library(FastLZeroSpikeInference)

idx = readRDS("../Data/Time_windows/indices.RDS")

idx_to_run = c(4,17,36,77)

n_w = 0
n_window = idx_to_run[n_w+1]

for(n_window in idx_to_run) {
  n_w = n_w+1
  cat(paste("Running window",n_w,"out of",length(idx_to_run)," - ID:",n_window,"\n"))
  
  filename = paste0("../Data/Time_windows/calcium_window", n_window, ".RDS")
  calcium = readRDS(filename)
  
  TT = nrow(calcium)
  N = ncol(calcium)

  estimated_calcium = matrix(0, N, TT)
  estimated_spikes = matrix(0, N, TT)
  estimated_amplitudes = matrix(0,N,TT)
  
  dd = apply(calcium, 2, diff)
  sigma2_start = var(calcium[calcium <  quantile(dd, .85)])
  
  #-----------------------# #-----------------------# 
  #      apply JW model for spike detection         #
  #-----------------------# #-----------------------# 
  lambdas = seq(0.01,10, length.out = 500)
  gam = 0.96
  sd_trace = sqrt(sigma2_start)
  
  # select the parameter lambda so that the number of detected spikes smaller than 1 estimated sd is zero
  
  for(j in 1:N){
    fit = sapply(lambdas,
                 function(lam) list(FastLZeroSpikeInference::estimate_spikes(dat = calcium[,j],
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
  
  
    fitL0 = FastLZeroSpikeInference::estimate_spikes(dat = calcium[,j],
                                                     gam = gam, lambda = lambda,
                                                     estimate_calcium = T)
    estimated_spikes[j, fitL0$spikes] = 1 
    estimated_calcium[j,] = fitL0$estimated_calcium
    
    for(t in which(estimated_spikes[j,]==1)){
      estimated_amplitudes[j,t] = estimated_calcium[j,t]-estimated_calcium[j,t-1]
      if(estimated_amplitudes[j,t] <0) { estimated_amplitudes[j,t] = 0 }
    }
  
  }
  

  name = paste0("02_data_analysis/04_comparison_two-step/results/estimated_spikes_JW_window", n_window, ".RDS")
  saveRDS(estimated_spikes, file = name)
  
  name = paste0("02_data_analysis/04_comparison_two-step/results/estimated_calcium_JW_window", n_window, ".RDS")
  saveRDS(estimated_calcium, file = name)
  
  name = paste0("02_data_analysis/04_comparison_two-step/results/estimated_amplitudes_JW_window", n_window, ".RDS")
  saveRDS(estimated_amplitudes, file = name)
  
  rm(estimated_spikes)
  rm(estimated_calcium)
  rm(fit)
  rm(fitL0)
  gc()

}
