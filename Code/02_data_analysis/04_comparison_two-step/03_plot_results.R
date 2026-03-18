#----------# #----------# #----------# #----------# #----------# #----------# #----------#
#----------#   THIS SCRIPT EXTRACTS ESTIMATED QUANTITIES AND PLOTS RESULTS    #----------#
#----------# #----------# #----------# #----------# #----------# #----------# #----------#


library(ggplot2)
library(ggpubr)
library(mclust)




#--------------------# #--------------------# 
##          Auxiliary functions            ##
#--------------------# #--------------------# 

sort_labels_by_size = function(est_cluster) {
  est_cluster = est_cluster + 100
  idc = sort(table(est_cluster), decreasing=T)
  for(i in 1:length(unique(est_cluster))){
    est_cluster[ est_cluster== as.numeric(attr(idc[i], "names")) ] = i
  }
  return(est_cluster)
}



#--------------------# #--------------------# 
##             Import data                 ##
#--------------------# #--------------------# 

loc_neurons <- readRDS("../Data/M3424F_loc_neurons.RDS")
indices_win = c(4,17,36,77)

indind <- 0


for(n_window in indices_win)
{
  indind <- indind+1
  cat(paste("Run", indind, "of", length(indices_win),"\n"))
  
  #----------# import run  
  name = paste0("02_data_analysis/04_comparison_two-step/results/estimated_spikes_JW_window", n_window, ".RDS")
  estimated_spikes = readRDS(file = name)
  
  name = paste0("02_data_analysis/04_comparison_two-step/results/estimated_calcium_JW_window", n_window, ".RDS")
  estimated_calcium = readRDS(file = name)
  
  name = paste0("02_data_analysis/04_comparison_two-step/results/estimates_clusterkmeans_window", n_window, ".RDS")
  estimated_cluster_kmeans = readRDS(file = name)
  
  #----------# import data
  cat("import data\n")
  filename = paste0("../Data/Time_windows/calcium_window", n_window, ".RDS")
  calcium = readRDS(filename)
  
  
  TT = nrow(calcium)
  N = ncol(calcium)
  
  estimated_cluster_kmeans = sort_labels_by_size(estimated_cluster_kmeans)
  table(estimated_cluster_kmeans)

  if(n_window == 17){
    # rename the clusters to improve visualization
    estimated_cluster_kmeans[estimated_cluster_kmeans==3] = 103
    estimated_cluster_kmeans[estimated_cluster_kmeans==4] = 3
    estimated_cluster_kmeans[estimated_cluster_kmeans==103] = 4
  }

  estimated_amplitudes = matrix(0,N,TT)
  for(i in 1:N){
    for(t in which(estimated_spikes[i,]==1)){
      estimated_amplitudes[i,t] = estimated_calcium[i,t]-estimated_calcium[i,t-1]
    }
  }
  estimated_amplitudes[estimated_amplitudes<0]=0
  
  par(mfrow=c(3,2))
  for(j in 1:6){
  plot(calcium[,j], type="l", ylim=c(-2,max(calcium[,j])+2))
  points(which(estimated_amplitudes[j,]>0), estimated_amplitudes[j,][which(estimated_amplitudes[j,]>0)] )
  lines(estimated_calcium[j,],col=2)
  }
  
  tapply(rowMeans(estimated_spikes), estimated_cluster_kmeans, mean)
  
  #---------# #-----------# #---------# #---------# #-----------# #---------# 
  #---------#    PLOT HEATMAP OF ACTIVATION SORTED BY CLUSTER     #---------#
  #---------# #-----------# #---------# #---------# #-----------# #---------# 
  
  cat("heatmap of the activations sorted by cluster\n")
  first_singleton = as.numeric(attr(table(estimated_cluster_kmeans), "names")[table(estimated_cluster_kmeans)==1][1]) 
  
  #####
  if(is.na(first_singleton)){
    first_singleton <- 999
  }
  
  
  order_series = sort(estimated_cluster_kmeans[(estimated_cluster_kmeans>1)&(estimated_cluster_kmeans<first_singleton)], index.return=T)
  annotations = c(0,which(diff(order_series$x)>0)) + c(diff( c(0,which(diff(order_series$x)>0)) )/2, 1) + 0.5
  
  mAA <- estimated_amplitudes[(estimated_cluster_kmeans>1)&(estimated_cluster_kmeans<first_singleton),]
  mAA <- mAA[order_series$ix,]
  meltM =reshape2::melt(t(mAA))
  meltM$time = meltM$Var1/15
  
  heatmap_spike = 
    ggplot() + 
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
  
  heatmap_spike
  
  filename = paste0("02_data_analysis/04_comparison_two-step/output_images/two-step_heatmap_spikes_win", n_window, ".pdf")
  ggsave(filename, heatmap_spike, width = 8, height = 6)
  filename = paste0("02_data_analysis/04_comparison_two-step/output_images/two-step_heatmap_spikes_win", n_window, ".png")
  ggsave(filename, heatmap_spike, width = 8, height = 6)

  
  
  #---------# #-----------# #---------# #---------# #-----------# #---------# 
  #---------#       PLOT CALCIUM TRACES SORTED BY CLUSTER         #---------#
  #---------# #-----------# #---------# #---------# #-----------# #---------# 
  
  idx = sort(c(estimated_cluster_kmeans), index.return = T)$ix
  idx = idx[ (sort(c(estimated_cluster_kmeans))>1) & 
               (sort(c(estimated_cluster_kmeans))<first_singleton) ]
  
  length(estimated_cluster_kmeans[idx]) # how many active neurons to plot?
  length(estimated_cluster_kmeans[idx])/8*3  # how many neurons to plot in a full column
  
  
  # rescale the series 
  calcium_active = calcium[,idx]
  calcium_active = apply(calcium_active, 2,function(x) x/25)
  
  TT = nrow(calcium_active)
  n = ncol(calcium_active)
  
  calcium_active = calcium_active + matrix(rep(1:n, TT), TT, n, byrow=T)
  calcium_active = reshape2::melt((calcium_active))
  calcium_active$time = rep(1:TT, n)/15
  calcium_active$neuron = as.factor(as.numeric(calcium_active$Var2))
  calcium_active$spike = as.numeric(reshape2::melt((estimated_spikes[idx,]))$value)
  calcium_active$Cluster = as.factor(rep(estimated_cluster_kmeans[idx], each = TT))
  
  calcium_active = rev(calcium_active)
  
  colors = c("1" = "black", "2" = "#9E0142", "3"="gold3", "4"="#009eb0", 
             "5"="#0045bd",
             "6"= "#5ba300",
             "7"="chocolate1", 
             "8"="#ef2578", "9"="#532257", "10"="#feb991", "11" = "grey", "12" = "grey27"
  )
  
  cond1 = as.numeric(calcium_active$neuron) <= ceiling(length(estimated_cluster_kmeans[idx])/8*3)
  
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
  
  
  
  cond2 = (as.numeric(calcium_active$neuron) > ceiling(length(estimated_cluster_kmeans[idx])/8*3) )& 
    (as.numeric(calcium_active$neuron) < ceiling(length(estimated_cluster_kmeans[idx])/8*3)*2 )
  
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
  
  
  cond3 = as.numeric(calcium_active$neuron) >= ceiling(length(estimated_cluster_kmeans[idx])/8*3)*2
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
  
  
  dataf = data.frame(x = loc_neurons[,1], y = loc_neurons[,2], Cluster = factor(estimated_cluster_kmeans))
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
  g
  
  
  filename = paste0("02_data_analysis/04_comparison_two-step/output_images/two-step_clustered_series_locations_win", n_window, ".pdf")
  ggsave(filename, g, width = 8, height = 9)
  filename = paste0("02_data_analysis/04_comparison_two-step/output_images/two-step_clustered_series_locations_win", n_window, ".png")
  ggsave(filename, g, width = 8, height = 9)
  
}
