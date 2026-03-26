
#---------# #-----------# #---------# #---------# #-----------# #---------# 
#---------#     COMPUTE AND PLOT WINDOW-SPECIFIC QUANTITIES     #---------#
#---------# #-----------# #---------# #---------# #-----------# #---------# 

# This script produces the plots contained in the output_images folder:
# - heatmap of the spike amplitudes (Figure 3 in the main paper)
# - calcium traces sorted and colored by cluster (Figure 4 in the main paper)
# - trajectories of the spike probabilities, i.e., the transformed Gaussian processes (Figure 5 in the main paper)

# To avoid downloading all the runs of the Gibbs sampler, you can set load_precomputed = TRUE.
# In this case, the script loads some pre-computed quantities:
#    - df_GP_win#.RDS
#    - est_cluster_neurons_win#.RDS
#    - estimated_spikes_win#.RDS
#    - meltM_win#.RDS
#
# If you prefer running it "from scratch", set the logical variable load_precomputed to FALSE.

# Also in this case, if you wish to reproduce the analyses only for the subset of windows 
# reported in the paper, you may set load_precomputed = FALSE and run_on_subset = TRUE.
# On the contrary, if you wish to run the code on all the time windows, set run_on_subset = FALSE.
# The script then executes a for loop over the selected time windows, and automatically
# saves the output plots into the folder as pdf and png images.


load_precomputed = TRUE
run_on_subset = TRUE


library(ggplot2)
library(salso)
library(TeachingDemos)
library(ggpubr)
library(mclust)


loc_neurons <- readRDS("../Data/M3424F_loc_neurons.RDS")
indices_win <- readRDS("../Data/Time_windows/indices.RDS")


#--------------------# #--------------------# 
##          Auxiliary functions            ##
#--------------------# #--------------------# 

# function that sorts the cluster by decreasing cluster size
sort_labels_by_size = function(est_cluster) {
  est_cluster = est_cluster + 100
  idc = sort(table(est_cluster), decreasing=T)
  for(i in 1:length(unique(est_cluster))){
    est_cluster[ est_cluster== as.numeric(attr(idc[i], "names")) ] = i
  }
  return(est_cluster)
}

# function to compute the false discovery rate
FDRk = function(k, PPS){
  logi_tmp = PPS>k
  num = sum((1-PPS) * logi_tmp)
  denom = sum(logi_tmp)
  return(num/denom)
}


#--------------------# #--------------------# #--------------------# #--------------------# 

if(run_on_subset) {
  idx_to_run = c(4,17,36,77)
} else {
  idx_to_run = indices_win
}


indind <- 0

for(n_window in idx_to_run)
{
  indind <- indind+1
  cat(paste("Run", indind, "of", length(idx_to_run),"\n"))
    
  #----------# import run if required
  if(!load_precomputed){
    cat("import run\n")
    filename = paste0("02_data_analysis/01_bSCDC_individual_trials/results/run_gibbs_window", n_window, "_alow_1.RDS")
    out = readRDS(filename)
  }
  
  #----------# import data
  cat("import data\n")
  filename = paste0("../Data/Time_windows/calcium_window", n_window, ".RDS")
  calcium = readRDS(filename)
  
  
  TT = nrow(calcium)
  N = ncol(calcium)
  
  
  #----------# #----------# #----------# #----------#
  #----------# ESTIMATE DETECTED SPIKES  #----------#
  #----------# #----------# #----------# #----------#
  
  cat("POSITIVE SPIKES\n")
  
  if(load_precomputed) {
    
    filename = paste0("02_data_analysis/01_bSCDC_individual_trials/output_RDS/estimated_spikes_win", n_window, ".RDS")
    estimated_spikes = readRDS(file=filename)
    
  } else {
    
    # PPS is a matrix where each cell is the spike probability
    PPS = matrix(0, TT, N)
    for(i in 1:N) {
      PPS[,i] = apply(out$AA[,i,1:length(out$gamma)], 1, function(x) mean(x>0))
    }
    
    # we need to select a threshold so that we identify a spike if P(spike)>threshold
    # we fix the false discovery rate to be below 0.05
    FDRk = Vectorize(FDRk, "k")
    fdr_threshold = 0.05
    threshold = 0.9
    threshold = uniroot(function(k) FDRk(k, PPS)-fdr_threshold, c(0.001, 0.9))$root
    
    estimated_spikes = (PPS>threshold)
    
    filename = paste0("02_data_analysis/01_bSCDC_individual_trials/output_RDS/estimated_spikes_win", n_window, ".RDS")
    saveRDS(estimated_spikes, file=filename)
    
    rm(PPS)
    
  }

  
  
  
  
  #----------# #----------# #----------# #----------#
  #----------#  COMPUTE NEURONS CLUSTER  #----------#
  #----------# #----------# #----------# #----------#
  
  cat("cluster of neurons\n")
  
  if(load_precomputed) {
    
    filename = paste0("02_data_analysis/01_bSCDC_individual_trials/output_RDS/est_cluster_neurons_win", n_window, ".RDS")
    est_cluster_neurons = readRDS(file=filename)
    
  } else {
    
    est_cluster_neurons = salso(t(out$cluster_signal+1), maxNClusters = 100)
    est_cluster_neurons = sort_labels_by_size(est_cluster_neurons)
    # table(est_cluster_neurons)
    
    filename = paste0("02_data_analysis/01_bSCDC_individual_trials/output_RDS/est_cluster_neurons_win", n_window, ".RDS")
    saveRDS(est_cluster_neurons, file=filename)
    
  }
  
  
  
  
  
  #---------# #-----------# #---------# #---------# #-----------# #---------# 
  #---------#    PLOT HEATMAP OF ACTIVATION SORTED BY CLUSTER     #---------#
  #---------# #-----------# #---------# #---------# #-----------# #---------# 
  
  cat("heatmap of the activations sorted by cluster\n")
  first_singleton = as.numeric(attr(table(est_cluster_neurons), "names")[table(est_cluster_neurons)==1][1]) 
  
  if(is.na(first_singleton)){
    first_singleton <- 999
  }
  
  if(load_precomputed){
    
    filename = paste0("02_data_analysis/01_bSCDC_individual_trials/output_RDS/meltM_win", n_window, ".RDS")
    meltM = readRDS(file=filename)
    order_series = sort(est_cluster_neurons[(est_cluster_neurons>1)&(est_cluster_neurons<first_singleton)], index.return=T)
    annotations = c(0,which(diff(order_series$x)>0)) + c(diff( c(0,which(diff(order_series$x)>0)) )/2, 1) + 0.5
    
  } else {
    
    mAA <- t(apply(out$AA,c(1,2),mean))
    mAA = mAA[(est_cluster_neurons>1)&(est_cluster_neurons<first_singleton), ]
    order_series = sort(est_cluster_neurons[(est_cluster_neurons>1)&(est_cluster_neurons<first_singleton)], index.return=T)
    annotations = c(0,which(diff(order_series$x)>0)) + c(diff( c(0,which(diff(order_series$x)>0)) )/2, 1) + 0.5
    
    mAA = mAA[order_series$ix,]
    mAA = mAA[,2:TT]
    meltM =reshape2::melt(t(mAA))
    meltM$time = meltM$Var1/15
    
    filename = paste0("02_data_analysis/01_bSCDC_individual_trials/output_RDS/meltM_win", n_window, ".RDS")
    saveRDS(meltM, file=filename)
  }

  heatmap_spike = ggplot() + 
    geom_tile(aes( x = time, y = Var2, fill = value),
              data = meltM) +
    scale_fill_viridis_c("Amplitude", option = "magma", direction=-1, )  +
    geom_hline(yintercept = which(diff(order_series$x)>0)+0.5, 
               linetype="dashed", 
               color = "gray", linewidth=0.5)+
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill='transparent', color=NA), #transparent panel bg
      plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
      panel.grid.major = element_blank(), #remove major gridlines
      panel.grid.minor = element_blank(), #remove minor gridlines
      panel.border = element_rect(color = "darkgray ", fill=NA),
      axis.line.x.bottom = element_line(color="gray"),
      legend.position = "right",
      legend.background = element_rect(fill='transparent', color = 'transparent'), #transparent legend bg
      legend.box.background = element_rect(fill='transparent', color = 'transparent'), #transparent legend panel
      legend.text = element_text(size=10),
      strip.text = element_text(size=10),
      strip.background = element_rect( fill=NA, color="gray76" ),
      text = element_text(size=12)
    )+
    xlab("Time (seconds)")  +
    ylab("Neurons") +
    annotate("text", x = -0.2, y = annotations, label = as.character(2:(length(annotations)+1))) + scale_y_reverse() 
  
  # heatmap_spike
  
  filename = paste0("02_data_analysis/01_bSCDC_individual_trials/output_images/heatmap_spikes_win", n_window, ".pdf")
  ggsave(filename, heatmap_spike, width = 8, height = 6)
  filename = paste0("02_data_analysis/01_bSCDC_individual_trials/output_images/heatmap_spikes_win", n_window, ".png")
  ggsave(filename, heatmap_spike, width = 8, height = 6)
  
  
  
  
  #---------# #-----------# #---------# #---------# #-----------# #---------# 
  #---------#  PLOT CLUSTER-SPECIFIC LATENT SPIKE PROBABILITIES   #---------#
  #---------# #-----------# #---------# #---------# #-----------# #---------# 
  
  cat("estimates latent Gaussian processes\n")
  
  if(load_precomputed) {
    
    filename = paste0("02_data_analysis/01_bSCDC_individual_trials/output_RDS/df_GP_win", n_window, ".RDS")
    df_GP = readRDS(file=filename)
    
  } else {
    
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
    ind_patch <- min(min(7,first_singleton-1),max(est_cluster_neurons))
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
  }
  
  plot_GPs = ggplot(data = df_GP, aes(x = time, y = y, color=clus)) + 
    geom_ribbon(aes(ymin = lower, ymax = upper),
                alpha=0.3, col = "turquoise4", fill = "turquoise4", lwd= 0.3) +
    geom_line(aes(y = y), col = "turquoise4", lwd = 1 ) +
    theme_minimal() +
    theme(
      panel.spacing = unit(1, "lines"),
      panel.background = element_rect(fill='transparent', color=NA), #transparent panel bg
      plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
      # panel.grid.major = element_blank(), #remove major gridlines
      panel.grid.minor = element_blank(), #remove minor gridlines
      panel.border = element_rect(color = "darkgray ", fill=NA),
      # axis.line.y.left = element_line(color="gray"),
      axis.line.x.bottom = element_line(color="gray"),
      #
      legend.position = "bottom",
      legend.background = element_rect(fill='transparent', color = 'white'), #transparent legend bg
      legend.box.background = element_rect(fill='transparent'), #transparent legend panel
      legend.text = element_text(size=10),
      strip.text = element_text(size=10),
      strip.background = element_rect( fill=NA, color="gray" )
    )+
    # scale_y_continuous(breaks = c(0.0, 0.1, 0.2,0.3),  limits = c(0,0.3)) +
    ylab("Spike probability")  + xlab("Time (seconds)") + 
    facet_wrap((label_clus)~., ncol=2)
  # plot_GPs
  
  filename = paste0("02_data_analysis/01_bSCDC_individual_trials/output_images/plot_GPs_win", n_window, ".pdf")
  ggsave(filename, plot_GPs, width = 8, height = 4.5)
  filename = paste0("02_data_analysis/01_bSCDC_individual_trials/output_images/plot_GPs_win", n_window, ".png")
  ggsave(filename, plot_GPs, width = 8, height = 4.5)
  
  
  rm(df_GP)
  rm(plot_GPs)
  
  
  
  
  #---------# #-----------# #---------# #---------# #-----------# #---------# 
  #---------#       PLOT CALCIUM TRACES SORTED BY CLUSTER         #---------#
  #---------# #-----------# #---------# #---------# #-----------# #---------# 
  
  idx = sort(c(est_cluster_neurons), index.return = T)$ix
  idx = idx[ (sort(c(est_cluster_neurons))>1) & 
             (sort(c(est_cluster_neurons))<first_singleton) ]
  
  length(est_cluster_neurons[idx]) # how many active neurons to plot?
  length(est_cluster_neurons[idx])/8*3  # how many neurons to plot in a full column
  
  
  # rescale the series 
  calcium_active = calcium[,idx]
  calcium_active = apply(calcium_active, 2,function(x) x/25)
  
  TT = nrow(calcium_active)
  n = ncol(calcium_active)
  
  calcium_active = calcium_active + matrix(rep(1:n, TT), TT, n, byrow=T)
  calcium_active = reshape2::melt((calcium_active))
  calcium_active$time = rep(1:TT, n)/15
  calcium_active$neuron = as.factor(as.numeric(calcium_active$Var2))
  calcium_active$spike = as.numeric(reshape2::melt((estimated_spikes[,idx]))$value)
  calcium_active$Cluster = as.factor(rep(est_cluster_neurons[idx], each = TT))
  
  calcium_active = rev(calcium_active)
  
  colors = c("1" = "black", "2" = "#9E0142", "3"="gold3", "4"="#009eb0", "5"="#0045bd","6"= "#5ba300",
             "7"="chocolate1", "8"="#ef2578", "9"="#532257", "10"="#feb991", "11" = "grey", "12" = "grey27"
  )
  
  cond1 = as.numeric(calcium_active$neuron) <= ceiling(length(est_cluster_neurons[idx])/8*3)
  
  length(unique(calcium_active$neuron[cond1]))
  subset_active = calcium_active[cond1,][calcium_active$spike[cond1]==1,]
  
  data_sub = calcium_active[cond1,]
  data_sub$rev_neuron = rev(data_sub$neuron)
  
  p1 = ggplot(data = data_sub, aes(x = time, y = value, color = Cluster)) +
    geom_ribbon(aes(x = time, ymax = value, ymin =as.numeric(neuron), fill=Cluster, color=Cluster, group = rev_neuron), alpha=0.08, lwd=0.1,
                show.legend=TRUE)+
    geom_line(aes(x = time, y = value, color=Cluster, group = rev_neuron), alpha=1, lwd=0.3,
              show.legend=TRUE)  +
    geom_point(data = subset_active, aes(x = time, y = as.numeric(neuron)), alpha = 0.35, size = 0.9,
               show.legend=TRUE) +
    scale_color_manual(values = colors, drop=FALSE,
                       breaks = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11","12")) +
    scale_fill_manual(values = colors,drop=FALSE,
                      breaks = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11","12")
    ) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill='transparent', color=NA), #transparent panel bg
      plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
      # panel.grid.major = element_blank(), #remove major gridlines
      panel.grid.minor = element_blank(), #remove minor gridlines
      panel.border = element_rect(color = "darkgray ", fill=NA),
      # axis.line.y.left = element_line(color="gray"),
      axis.line.x.bottom = element_line(color="gray"),
      #
      legend.position = "bottom",
      legend.background = element_rect(fill='transparent', color = 'transparent'), #transparent legend bg
      legend.box.background = element_rect(fill='transparent', color = 'transparent'), #transparent legend panel
      legend.text = element_text(size=10),
      strip.text = element_text(size=10),
      strip.background = element_rect( fill=NA, color="gray" )
    )+
    guides(colour = guide_legend(nrow = 1)) +
    xlab("Time (seconds)")  +
    ylab("Neuron") 
  
  
  
  cond2 = (as.numeric(calcium_active$neuron) > ceiling(length(est_cluster_neurons[idx])/8*3) )& 
    (as.numeric(calcium_active$neuron) < ceiling(length(est_cluster_neurons[idx])/8*3)*2 )
  
  length(unique(calcium_active$neuron[cond2]))
  
  subset_active = calcium_active[cond2,][calcium_active$spike[cond2]==1,]
  data_sub = calcium_active[cond2,]
  
  data_sub$rev_neuron = rev(data_sub$neuron)
  
  p2 = ggplot(data = data_sub, aes(x = time, y = value, color = Cluster)) +
    geom_ribbon(aes(x = time, ymax = value, ymin =as.numeric(neuron), fill=Cluster, color=Cluster, group = rev_neuron), alpha=0.08, lwd=0.1,
                show.legend=TRUE)+
    geom_line(aes(x = time, y = value, color=Cluster, group = rev_neuron), alpha=1, lwd=0.3,
              show.legend=TRUE)  +
    geom_point(data = subset_active, aes(x = time, y = as.numeric(neuron)), alpha = 0.35, size = 0.9,
               show.legend=TRUE) +
    scale_color_manual(values = colors,
                       breaks = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11","12"),
                       drop = FALSE ) +
    scale_fill_manual(values = colors,
                      breaks = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11","12"),
                      drop = FALSE ) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill='transparent', color=NA), #transparent panel bg
      plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
      # panel.grid.major = element_blank(), #remove major gridlines
      panel.grid.minor = element_blank(), #remove minor gridlines
      panel.border = element_rect(color = "darkgray ", fill=NA),
      # axis.line.y.left = element_line(color="gray"),
      axis.line.x.bottom = element_line(color="gray"),
      #
      legend.position = "bottom",
      legend.background = element_rect(fill='transparent', color = 'white'), #transparent legend bg
      legend.box.background = element_rect(fill='transparent'), #transparent legend panel
      legend.text = element_text(size=10),
      strip.text = element_text(size=10),
      strip.background = element_rect( fill=NA, color="gray" )
    )+
    guides(colour = guide_legend(nrow = 1)) +
    xlab("Time (seconds)")  +
    ylab("")
  
  
  cond3 = as.numeric(calcium_active$neuron) >= ceiling(length(est_cluster_neurons[idx])/8*3)*2
  length(unique(calcium_active$neuron[cond3]))
  subset_active = calcium_active[cond3,][calcium_active$spike[cond3]==1,]
  data_sub = calcium_active[cond3,]
  
  data_sub$rev_neuron = rev(data_sub$neuron)
  
  p3 = ggplot(data = data_sub, aes(x = time, y = value, color = Cluster)) +
    geom_ribbon(aes(x = time, ymax = value, ymin =as.numeric(neuron), fill=Cluster, color=Cluster, group = rev_neuron), alpha=0.08, lwd=0.1,
                show.legend=TRUE)+
    geom_line(aes(x = time, y = value, color=Cluster, group = rev_neuron), alpha=1, lwd=0.3,
              show.legend=TRUE)  +
    geom_point(data = subset_active, aes(x = time, y = as.numeric(neuron)), alpha = 0.35, size = 0.9,
               show.legend=TRUE) +
    scale_color_manual(values = colors,
                       breaks = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11","12"),
                       drop = FALSE ) +
    scale_fill_manual(values = colors,
                      breaks = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11","12"),
                      drop = FALSE ) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill='transparent', color=NA), #transparent panel bg
      plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
      # panel.grid.major = element_blank(), #remove major gridlines
      panel.grid.minor = element_blank(), #remove minor gridlines
      panel.border = element_rect(color = "darkgray ", fill=NA),
      # axis.line.y.left = element_line(color="gray"),
      axis.line.x.bottom = element_line(color="gray"),
      #
      legend.position = "bottom",
      legend.background = element_rect(fill='transparent', color = 'white'), #transparent legend bg
      legend.box.background = element_rect(fill='transparent'), #transparent legend panel
      legend.text = element_text(size=10),
      strip.text = element_text(size=10),
      strip.background = element_rect( fill=NA, color="gray" )
    )+
    guides(colour = guide_legend(nrow = 1)) +
    xlab("Time (seconds)")  + 
    ylab("") 
  
  
  dataf = data.frame(x = loc_neurons[,1], y = loc_neurons[,2], Cluster = factor(est_cluster_neurons))
  p4 = ggplot(dataf, aes(x = x, y=y))+
    scale_fill_manual(values = colors,
                      breaks = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11","12"),
                      drop = FALSE ) +
    geom_point(aes(fill=Cluster),
               colour="black", size = 2.4, pch=21, alpha=0.8) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill='transparent', color=NA), #transparent panel bg
      plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
      # panel.grid.major = element_blank(), #remove major gridlines
      panel.grid.minor = element_blank(), #remove minor gridlines
      panel.border = element_rect(color = "darkgray ", fill=NA),
      # axis.line.y.left = element_line(color="gray"),
      axis.line.x.bottom = element_line(color="gray"),
      #
      legend.direction='vertical',
      legend.position = "right",
      legend.background = element_rect(fill='transparent', color = 'white'), #transparent legend bg
      legend.box.background = element_rect(fill='transparent'), #transparent legend panel
      legend.text = element_text(size=10),
      strip.text = element_text(size=10),
      strip.background = element_rect( fill=NA, color="gray" )
    )+
    xlab("x coordinate")  +
    ylab("y coordinate")
  
  
  g <- ggarrange(p1, p2, 
                 ggarrange(p3,p4, nrow=2, ncol=1, legend = F, heights = c(1.9,1)), 
                 ncol=3, nrow=1,
                 common.legend = TRUE, legend="bottom")
  # g
  
  filename = paste0("02_data_analysis/01_bSCDC_individual_trials/output_images/clustered_series_locations_win", n_window, ".pdf")
  ggsave(filename, g, width = 8, height = 9)
  filename = paste0("02_data_analysis/01_bSCDC_individual_trials/output_images/clustered_series_locations_win", n_window, ".png")
  ggsave(filename, g, width = 8, height = 9)
  
  if(!load_precomputed){rm(out)}
  rm(calcium, calcium_active)
  rm(data_sub, dataf)
  rm(estimated_spikes, est_cluster_neurons)
  rm(meltM)
  rm(p1,p2,p3,p4)
  rm(cond1,cond2,cond3)
  rm(subset_active)
  rm(g, heatmap_spike, order_series)
  gc()
}
