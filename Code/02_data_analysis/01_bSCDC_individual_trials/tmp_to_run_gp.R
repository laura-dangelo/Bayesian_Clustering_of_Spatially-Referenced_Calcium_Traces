

library(ggplot2)
library(salso)
library(TeachingDemos)
library(ggpubr)
library(mclust)


loc_neurons <- readRDS("../Data/M3424F_loc_neurons.RDS")
indices_win <- readRDS("../Data/Time_windows/indices.RDS")



#--------------------# #--------------------# #--------------------# #--------------------# 

idx_to_run = indices_win

indind <- 0

for(n_window in idx_to_run)
{
  indind <- indind+1
  cat(paste("Run", indind, "of", length(idx_to_run),"\n"))
  
  #----------# import run if required
  filename = paste0("02_data_analysis/01_bSCDC_individual_trials/results/run_gibbs_window", n_window, "_alow_1.RDS")
  out = readRDS(filename)
  
  
  #----------# import data
  cat("import data\n")
  filename = paste0("../Data/Time_windows/calcium_window", n_window, ".RDS")
  calcium = readRDS(filename)
  
  
  TT = nrow(calcium)
  N = ncol(calcium)
  
  filename = paste0("02_data_analysis/01_bSCDC_individual_trials/output_RDS/est_cluster_neurons_win", n_window, ".RDS")
  est_cluster_neurons = readRDS(file=filename)
  

  # select the mcmc iterations with a cluster similar to the estimate
  selected_rows = c()
  j=1
  ARIS = rev(seq(0.3,1, by=0.1))
  while((length(selected_rows)<500)&j<length(ARIS)){
    id = sapply(1:length(out$gamma), function(x) adjustedRandIndex(est_cluster_neurons, (out$cluster_signal)[,x] ) > ARIS[j] )
    selected_rows = which(id==T)
    j=j+1
  }
  
  mcmc_GP = out$latent_signal[,,selected_rows]
  ind_patch <- min(first_singleton-1,max(est_cluster_neurons))
  clus = 2:ind_patch
  df_GP = data.frame("clus"=sort(rep(clus,TT)), "time"=rep(1:TT,length(clus)), 
                       "y"=0, "lower" = 0, "upper" = 0)
    
  for(i in 1:length(clus)){
    obs_cl = which(est_cluster_neurons==clus[i])[1]
    seq_cl = out$cluster_signal[obs_cl,selected_rows]+1
      
    est_GP = matrix(0, length(selected_rows), TT)
    for(iter in 1:length(selected_rows)) est_GP[iter,] = 
      pnorm( out$latent_signal[, seq_cl[iter], selected_rows[iter]] )
    hpdGP = apply(est_GP, 2, function(x) emp.hpd(x))
      
    df_GP[(df_GP$clus == clus[i]),]$y = colMeans(est_GP)
    df_GP[df_GP$clus == clus[i],]$lower = hpdGP[1,]
    df_GP[df_GP$clus == clus[i],]$upper = hpdGP[2,]
  }
    
  df_GP$label_clus = paste0("Cluster ", df_GP$clus)
  df_GP$time = df_GP$time/15
    
  filename = paste0("02_data_analysis/01_bSCDC_individual_trials/output_RDS/df_GP_win", n_window, ".RDS")
  saveRDS(df_GP, file=filename)
  rm(id,i,j,selected_rows, seq_cl, est_GP, mcmc_GP,hpdGP)

  rm(out)
  gc()
}
