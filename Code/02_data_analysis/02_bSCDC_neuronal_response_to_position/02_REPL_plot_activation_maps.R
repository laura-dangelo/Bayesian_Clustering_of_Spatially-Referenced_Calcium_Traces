
#---------# #-----------# #---------# #---------# #-----------# #---------# 
#---------#        PLOT THE MAP OF THE NEURONS' ACTIVITY        #---------#
#---------# #-----------# #---------# #---------# #-----------# #---------# 
# This script produces the plots of Figure 7 of the main paper, which show how the activity of individual neurons
# is affected by the mouse's position in the arena.

library(ggplot2)
library(tidyverse)
library(viridisLite)
library(ggpubr)
library(patchwork)

source("02_data_analysis/02_bSCDC_neuronal_response_to_position/00_auxiliary_functions_DONT_RUN.R")

data <- readRDS("../Data/data_binary_position.RDS")
WIND <- readRDS("02_data_analysis/02_bSCDC_neuronal_response_to_position/output_RDS/Windows_list.RDS")
idx <- readRDS("../Data/Time_windows/indices.RDS")
loc_neurons <- readRDS("../Data/M3424F_loc_neurons.RDS")



#---------# #---------# #---------# #---------# #---------# 
#   Pick neurons with a large co-clustering probability   #
#---------# #---------# #---------# #---------# #---------# 

est_cluster <- NULL
for(i in 1:length(idx)){
  n_window <- idx[i]
  filename <- paste0("02_data_analysis/01_bSCDC_individual_trials/output_RDS/est_cluster_neurons_win", n_window, ".RDS")
  est_cl <- readRDS(file=filename)
  est_cluster <-rbind(est_cluster, est_cl)
}
rm(est_cl)
dim(est_cluster)

cocluster_neurons <- matrix(1,ncol(est_cluster),ncol(est_cluster))
for(i in 2:nrow(cocluster_neurons)){
  for(j in 1:(i-1)) {
    cocluster_neurons[i,j] = cocluster_neurons[j,i] = sum(apply(est_cluster[,c(i,j)], 1, function(x) x[1]==x[2]))
  }
}

summary(c(cocluster_neurons)/80)
which(cocluster_neurons/80 > .8, arr.ind = TRUE)


neu1 = 312
neu2 = 288
neu3 = 232
neu4 = 118

pclust <- cocluster_neurons/80
pclust[neu1,neu2]
pclust[neu3,neu4]



#---------# #---------# #---------# #---------# 
#       Plot neurons' spiking activity        #
#---------# #---------# #---------# #---------# 

data$pos_binary[7867:(7867+256)] = 3
data$pos_binary[(7867+256):8378] = 4
idcirc = (data$pos1^2 + data$pos2^2)<1.05
win = unlist(apply(
  cbind(
    1:length(which(diff(data$pos_binary)!=0)),
    c(0,(which(diff(data$pos_binary)!=0)))[1:138],
    (which(diff(data$pos_binary)!=0))), 1, function(x) rep(x[1], x[3]-x[2])))
win = c(win, rep(139, 10870-10862))
data$win = win
data = data[,c(1:5,331,6:330)]


G11 <- plot_neuron_smooth_spatial(neu = neu1, data_df = data,
                                  WIND_list = WIND, idx_vec = idx,
                                  upper = .50,smooth = 5,npix = 150)
G12 <- plot_neuron_smooth_spatial(neu = neu2, data,WIND_list = WIND,idx_vec = idx,
                                  upper = .50,smooth = 5,npix = 150)
G11+G12

G21 <- plot_neuron_smooth_spatial(neu = neu3,data,WIND_list = WIND,idx_vec = idx,
                                  upper = .50,smooth = 5,npix = 150)
G22 <- plot_neuron_smooth_spatial(neu = neu4,data,WIND_list = WIND,idx_vec = idx,
                                  upper = .50,smooth = 5,npix = 150)
G21+G22

Q11 <- G11+xlab(" \n ")+ggtitle("Neuron A")
Q12 <- G12+ylab(" \n ")+xlab(" \n ")+ggtitle("Neuron B")
Q21 <- G21+ggtitle("Neuron C")
Q22 <- G22+ylab(" \n ")+ggtitle("Neuron D")



#---------# #---------# #---------# #---------# 
#  Plot neurons' locations in the hippocampus #
#---------# #---------# #---------# #---------# 

G13 <- plot_neuron_locations(neu1, neu2, c(23,24))
G23 <- plot_neuron_locations(neu3, neu4, c(22,25), cols = c("forestgreen","maroon2"))
Q13 <- G13 +ggtitle("Neurons A and B")+xlab(" \n ")
Q23 <- G23 +ggtitle("Neurons C and D")


G = ggarrange(Q11, Q12, Q13, Q21, Q22, Q23, widths = rep(1,6),
              ncol=3, nrow=2,
              common.legend = TRUE, legend="bottom")
G
ggsave(filename = "02_data_analysis/02_bSCDC_neuronal_response_to_position/output_images/neuronal_activation_smooth.pdf",
       width = 14, height = 14/3*2)
ggsave(filename = "02_data_analysis/02_bSCDC_neuronal_response_to_position/output_images/neuronal_activation_smooth.png",
       width = 14, height = 14/3*2)

