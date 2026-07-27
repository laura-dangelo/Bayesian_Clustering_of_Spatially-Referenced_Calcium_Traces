setwd("02_data_analysis/02_bSCDC_neuronal_response_to_position/")

source("00_auxiliary_functions_DONT_RUN.R")
data <- readRDS("../../../Data/data_binary_position.RDS")
WIND <- readRDS("../../02_data_analysis/02_bSCDC_neuronal_response_to_position/output_RDS/Windows_list.RDS")
idx <- readRDS("../../../Data/Time_windows/indices.RDS")
loc_neurons <- readRDS("../../../Data/M3424F_loc_neurons.RDS")
total_neurons_available <- nrow(loc_neurons)

est_cluster <- NULL
for(i in 1:length(idx)){
  n_window <- idx[i]
  filename <- paste0("../01_bSCDC_individual_trials/output_RDS/est_cluster_neurons_win", n_window, ".RDS")
  est_cl <- readRDS(file=filename)
  est_cluster <-rbind(est_cluster, est_cl)
}
rm(est_cl)

cocluster_neurons <- matrix(1,ncol(est_cluster),ncol(est_cluster))
for(i in 2:nrow(cocluster_neurons)){
  for(j in 1:(i-1)) {
    cocluster_neurons[i,j] = cocluster_neurons[j,i] = sum(apply(est_cluster[,c(i,j)], 1, function(x) x[1]==x[2]))
  }
}

# -------------------------------------------------------------------------


# add sequential window number to the calcium data
str(data)
data$pos_binary[7867:(7867+256)] = 3
data$pos_binary[(7867+256):8378] = 4
idcirc = (data$pos1^2 + data$pos2^2)<1.05
win = unlist(apply(
  cbind(
    1:length(which(diff(data$pos_binary)!=0)),
    c(0,(which(diff(data$pos_binary)!=0)))[1:138],
    (which(diff(data$pos_binary)!=0))), 1, function(x) rep(x[1], x[3]-x[2])))
win = c(win, rep(139, nrow(data)-length(win)))
data$win = win
data = data[,c(1:5,331,6:330)]


for(ii in 1:total_neurons_available){
V <- plot_neuron_smooth_spatial(neu = ii,
                           data_df = data,
                           WIND_list = WIND,
                           idx_vec = idx,
                           npix = 250,
                           upper = 1,smooth = 5)
ggsave(filename = paste0("output_images/all_maps_neurons_interpolation/interpolation_neuron", ii, ".pdf"),
       plot = V,
       device = cairo_pdf,
       width = 10, height = 10)
cat(paste("Done with neuron",ii,"out of", total_neurons_available,"\n"))
}
