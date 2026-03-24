
#---------# #-----------# #---------# #---------# #-----------# #---------# #---------# 
#---------#         COMPUTE QUATITIES FOR SENSITIVITY ANALYSIS A_LOW        #---------#
#---------# #-----------# #---------# #---------# #-----------# #---------# #---------# 

# First part of the sensitivity analysis reported in Section S2.1 of the Supplementary Material. 
# This script extract the inferred spike trains and cluster of neurons.

#  _________________________________
#  YOU CAN AVOID RUNNING THIS SCRIPT 
#  ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
# This script produces the RDS files available in the output_RDS folder; specifically:
#   - est_cluster_neurons_win17_alow_#.RDS
#   - estimated_spikes_win17_alow_#.RDS
# These quantities are then used in the next script, which compares the results for different values of a_low.

library(ggplot2)
library(salso)
library(TeachingDemos)
library(ggpubr)
library(mclust)

loc_neurons = readRDS("../Data/M3424F_loc_neurons.RDS")
indices_win = readRDS("../Data/Time_windows/indices.RDS")

# #----------#  pick a window
n_window = 17
a_low_seq = c(0, 0.5)

for(a_low in a_low_seq) {
  
  #----------#  import run
  filename = paste0("02_data_analysis/03_bSCDC_sensitivity_study/results/run_gibbs_window", n_window,"_alow_", a_low,".RDS")
  out = readRDS(filename)
  
  #----------#  import data
  filename = paste0("../Data/Time_windows/calcium_window", n_window, ".RDS")
  calcium = readRDS(filename)
  
  TT = nrow(calcium)
  N = ncol(calcium)
  

  #----------# #----------# #----------# #----------#
  #----------# ESTIMATE DETECTED SPIKES  #----------#
  #----------# #----------# #----------# #----------#
  
  # PPS is a matrix where each cell is the spike probability
  PPS = matrix(0, TT, N)
  for(i in 1:N) {
    PPS[,i] = apply(out$AA[,i,1:length(out$gamma)], 1, function(x) mean(x>0))
  }
  
  # we need to select a threshold so that we identify a spike if P(spike)>threshold
  # we fix the false discovery rate to be below 0.05
  FDRk = function(k, PPS){
    logi_tmp = PPS>k
    num = sum((1-PPS) * logi_tmp)
    denom = sum(logi_tmp)
    return(num/denom)
  }
  
  FDRk = Vectorize(FDRk, "k")
  
  fdr_threshold = 0.05
  threshold = 0.9
  threshold = uniroot(function(k) FDRk(k, PPS)-fdr_threshold, c(0.001, 0.9))$root
  estimated_spikes = (PPS>threshold)
  
  filename = paste0("02_data_analysis/03_bSCDC_sensitivity_study/output_RDS/estimated_spikes_win", n_window,"_alow_", a_low,".RDS")
  saveRDS(estimated_spikes, file=filename)
  
  rm(PPS)
  
  
  
  
  #----------# #----------# #----------# #----------#
  #----------#  COMPUTE NEURONS CLUSTER  #----------#
  #----------# #----------# #----------# #----------#
  
  est_cluster_neurons = salso(t(out$cluster_signal+1), maxNClusters = 100)
  
  sort_labels_by_size = function(est_cluster) {
    est_cluster = est_cluster + 100
    idc = sort(table(est_cluster), decreasing=T)
    for(i in 1:length(unique(est_cluster))){
      est_cluster[ est_cluster== as.numeric(attr(idc[i], "names")) ] = i
    }
    return(est_cluster)
  }
  est_cluster_neurons = sort_labels_by_size(est_cluster_neurons)
  # table(est_cluster_neurons)
  
  filename = paste0("02_data_analysis/03_bSCDC_sensitivity_study/output_RDS/est_cluster_neurons_win", n_window,"_alow_", a_low,".RDS")
  saveRDS(est_cluster_neurons, file=filename)
  
  
  
  rm(calcium)
  rm(estimated_spikes, est_cluster_neurons)
  rm(out)
  gc()
}


filename = paste0("02_data_analysis/01_bSCDC_individual_trials/output_RDS/estimated_spikes_win", n_window, ".RDS")
estimated_spikes = readRDS(file=filename)
filename = paste0("02_data_analysis/03_bSCDC_sensitivity_study/output_RDS/estimated_spikes_win", n_window,"_alow_1.RDS")
saveRDS(estimated_spikes, file=filename)


filename = paste0("02_data_analysis/01_bSCDC_individual_trials/output_RDS/est_cluster_neurons_win", n_window, ".RDS")
est_cluster_neurons = readRDS(file=filename)
filename = paste0("02_data_analysis/03_bSCDC_sensitivity_study/output_RDS/est_cluster_neurons_win", n_window,"_alow_1.RDS")
saveRDS(est_cluster_neurons, file=filename)

