
#---------# #-----------# #---------# #---------# #-----------# #---------# #---------# 
#---------#                 RESULTS  SENSITIVITY ANALYSIS PSBP              #---------#
#---------# #-----------# #---------# #---------# #-----------# #---------# #---------# 

# Second part of the sensitivity analysis reported in Section S2.1 of the Supplementary Material. 
# This script produces the quantities of interest reported in the article. 

library(plyr) 
library(ggplot2)
library(mclust)

estimates_sigma2 = numeric(3)
estimates_tau2 = numeric(3)
estimates_gamma = numeric(3)
firing_rates = numeric(3)

# #----------#  pick a window
n_window = 17
sigma2_alpha_seq = c(0.5, 0.1, 0.01)



for(i in 1:length(sigma2_alpha_seq)) {
  #----------#  import run
  filenameout = paste0("02_data_analysis/03_bSCDC_sensitivity_study/results/run_gibbs_window", n_window,"_theta_", 1/sigma2_alpha_seq[i],".RDS")
  out = readRDS(file = filenameout)

  # posterior estimates of the parameters
  estimates_sigma2[i] = mean(out$sigma2)
  estimates_tau2[i] = mean(out$tau2)
  estimates_gamma[i] = mean(out$gamma)
  rm(out)
  
  # firing rates
  filename = paste0("02_data_analysis/03_bSCDC_sensitivity_study/output_RDS/estimated_spikes_win", n_window,"_theta_", 1/sigma2_alpha_seq[i],".RDS")
  estimated_spikes = readRDS(file=filename)
  firing_rates[i] = mean(estimated_spikes)
  
}

estimates_sigma2
estimates_tau2
estimates_gamma
firing_rates


filename = paste0("02_data_analysis/03_bSCDC_sensitivity_study/output_RDS/est_cluster_neurons_win", n_window,"_theta_", 1/sigma2_alpha_seq[1],".RDS")
est_cluster_neurons05 = readRDS(file=filename)

filename = paste0("02_data_analysis/03_bSCDC_sensitivity_study/output_RDS/est_cluster_neurons_win", n_window,"_theta_", 1/sigma2_alpha_seq[2],".RDS")
est_cluster_neurons01 = readRDS(file=filename)

filename = paste0("02_data_analysis/03_bSCDC_sensitivity_study/output_RDS/est_cluster_neurons_win", n_window,"_theta_", 1/sigma2_alpha_seq[3],".RDS")
est_cluster_neurons001 = readRDS(file=filename)


sort_labels_by_size = function(est_cluster) {
  est_cluster = est_cluster + 100
  idc = sort(table(est_cluster), decreasing=T)
  for(i in 1:length(unique(est_cluster))){
    est_cluster[ est_cluster== as.numeric(attr(idc[i], "names")) ] = i
  }
  return(est_cluster)
}
est_cluster_neurons05 = sort_labels_by_size(est_cluster_neurons05)
est_cluster_neurons01 = sort_labels_by_size(est_cluster_neurons01)
est_cluster_neurons001 = sort_labels_by_size(est_cluster_neurons001)



ARIs = matrix(1,3,3)

ARIs[1,2] = ARIs[2,1] = adjustedRandIndex(est_cluster_neurons05, est_cluster_neurons01)
ARIs[1,3] = ARIs[3,1] = adjustedRandIndex(est_cluster_neurons05, est_cluster_neurons001)

ARIs[2,3] = ARIs[3,2] = adjustedRandIndex(est_cluster_neurons01, est_cluster_neurons001)

colnames(ARIs) = 1/sigma2_alpha_seq
rownames(ARIs) = 1/sigma2_alpha_seq
ARIs


table(est_cluster_neurons05)
table(est_cluster_neurons01)
table(est_cluster_neurons001)


est_cluster_neurons01 = est_cluster_neurons01 + 100
est_cluster_neurons01[est_cluster_neurons01==102] = 4
est_cluster_neurons01[est_cluster_neurons01==103] = 3
est_cluster_neurons01[est_cluster_neurons01==104] = 2
est_cluster_neurons01[est_cluster_neurons01>100] = est_cluster_neurons01[est_cluster_neurons01>100]-100
table(est_cluster_neurons01)



est_cluster_neurons001 = est_cluster_neurons001 + 100
est_cluster_neurons001[est_cluster_neurons001==102] = 3
est_cluster_neurons001[est_cluster_neurons001==103] = 2
est_cluster_neurons001[est_cluster_neurons001==105] = 8
est_cluster_neurons001[est_cluster_neurons001==108] = 5
est_cluster_neurons001[est_cluster_neurons001>100] = est_cluster_neurons001[est_cluster_neurons001>100]-100
table(est_cluster_neurons001)



loc_neurons = readRDS("../Data/M3424F_loc_neurons.RDS")

colors = c("1" = "black", "2" = "#9E0142", "3"="gold3", "4"="#009eb0", "5"="#0045bd","6"= "#5ba300",
           "7"="chocolate1", "8"="#ef2578", "9"="#532257", "10"="#feb991", "11" = "grey", "12" = "grey27",
           "13"="grey56", "14"="lavender", "15"="blue3","16"= "forestgreen",
           "17"="darkslategrey", "18"="grey44"
)
dataf = data.frame(x = loc_neurons[,1], y = loc_neurons[,2], 
                   Cluster = factor(est_cluster_neurons001),
                   theta = rep("theta = 100", length(loc_neurons[,1]))
                   )
dataf = rbind(dataf, 
              data.frame(x = loc_neurons[,1], y = loc_neurons[,2], 
                         Cluster = factor(est_cluster_neurons01),
                         theta = rep("theta = 10", length(loc_neurons[,1]))
              ))
dataf = rbind(dataf, 
              data.frame(x = loc_neurons[,1], y = loc_neurons[,2], 
                         Cluster = factor(est_cluster_neurons05),
                         theta = rep("theta = 2", length(loc_neurons[,1]))
              ))

dataf$theta = factor(dataf$theta)

neworder = c("theta = 2","theta = 10","theta = 100")
 
dataf = arrange(transform(dataf,
                           theta=factor(theta,levels=neworder)),theta)


p4 = ggplot(dataf, aes(x = x, y=y))+
  scale_fill_manual(values = colors,
                    breaks = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11","12","13",
                               "14", "15", "16", "17", "18"),
                    drop = FALSE ) +
  geom_point(aes(fill=Cluster),
             colour="black", size = 3, pch=21, alpha=0.8) +
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
  ylab("y coordinate") + 
  facet_wrap(~theta)
p4

ggsave("02_data_analysis/03_bSCDC_sensitivity_study/output_images/cluster_PSB.pdf", p4, width = 9, height = 5)
ggsave("02_data_analysis/03_bSCDC_sensitivity_study/output_images/cluster_PSB.png", p4, width = 9, height = 5)
