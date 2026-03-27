
#---------# #-----------# #---------# #---------# #-----------# #---------# #---------# 
#---------#                 RESULTS  SENSITIVITY ANALYSIS A_LOW             #---------#
#---------# #-----------# #---------# #---------# #-----------# #---------# #---------# 

# First part of the sensitivity analysis reported in Section S2.1 of the Supplementary Material. 
# This script produces the quantities of interest reported in the article. 

library(mclust)

estimates_sigma2 = numeric(3)
estimates_tau2 = numeric(3)
estimates_gamma = numeric(3)
firing_rates = numeric(3)
ARIs = matrix(1,3,3)

# #----------#  pick a window
n_window = 17
a_low_seq = c(0, 0.5, 1)



for(i in 1:length(a_low_seq)) {

  filename = paste0("02_data_analysis/03_bSCDC_sensitivity_study/output_RDS/estimates_win", n_window,"_alow_", a_low_seq[i],".RDS")
  estimates = readRDS(file=filename)
    
  # posterior estimates of the parameters
  estimates_sigma2[i] = estimates$sigma2
  estimates_tau2[i] = estimates$tau2
  estimates_gamma[i] = estimates$gamma

    
  # firing rates
  filename = paste0("02_data_analysis/03_bSCDC_sensitivity_study/output_RDS/estimated_spikes_win", n_window,"_alow_", a_low_seq[i],".RDS")
  estimated_spikes = readRDS(file=filename)
  firing_rates[i] = mean(estimated_spikes)
  
}

estimates_sigma2
estimates_tau2
estimates_gamma
firing_rates


filename = paste0("02_data_analysis/03_bSCDC_sensitivity_study/output_RDS/est_cluster_neurons_win", n_window,"_alow_", 0,".RDS")
est_cluster_neurons0 = readRDS(file=filename)

filename = paste0("02_data_analysis/03_bSCDC_sensitivity_study/output_RDS/est_cluster_neurons_win", n_window,"_alow_", 0.5,".RDS")
est_cluster_neurons05 = readRDS(file=filename)

filename = paste0("02_data_analysis/03_bSCDC_sensitivity_study/output_RDS/est_cluster_neurons_win", n_window,"_alow_", 1,".RDS")
est_cluster_neurons1 = readRDS(file=filename)

ARIs[1,2] = ARIs[2,1] = adjustedRandIndex(est_cluster_neurons0, est_cluster_neurons05)
ARIs[1,3] = ARIs[3,1] = adjustedRandIndex(est_cluster_neurons0, est_cluster_neurons1)
ARIs[2,3] = ARIs[3,2] = adjustedRandIndex(est_cluster_neurons1, est_cluster_neurons05)
ARIs
