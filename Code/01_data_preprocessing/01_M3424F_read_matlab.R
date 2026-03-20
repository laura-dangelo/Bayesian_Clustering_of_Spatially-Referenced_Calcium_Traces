#-------# #-------# #-------# #-------# #-------# #-------# #-------# #-------# 
#
#        Import M3424F_data_togo_neuron_behav_multiTrials_072621.mat 
#        and extract individual quantities 
#
#-------# #-------# #-------# #-------# #-------# #-------# #-------# #-------# 

# This script imports the original .mat data, extracts the individual quantities of interest and
# saves them as separate RDS files.
# In particular, it saves to the Data folder the following files:
#   - M3424F_loc_neurons.RDS
#   - M3424F_time_calcium.RDS
#   - M3424F_time_mouse_position.RDS
#   - M3424F_calcium.RDS                                     
#   - M3424F_coord_mouse_position.RDS

#  _________________________________
#  YOU CAN AVOID RUNNING THIS SCRIPT 
#  ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
# In the Google Drive folder (https://drive.google.com/drive/folders/1-1xf57mZBc1usA-iCZGkp4KPF8oX_5mV?usp=sharing)
# the RDS files produced by this script are available in the "Data" folder.
# You can download these files and copy them in the Data folder of this repository.



library(rmatio)
data = read.mat(filename = "../../Data/M3424F_data_togo_neuron_behav_multiTrials_072621.mat")
str(data)


# list of 2:
# $behav
# $neuron_dat

index = 4
# str(data$behav[[1]][[index]]) # behavior mouse
# str(data$neuron_dat[[1]] ) # neurons activity



#-------# #-------# #-------# #-------# #-------# #-------# #-------# 
#-------# #-------#       Mouse behavior data     #-------# #-------# 
#-------# #-------# #-------# #-------# #-------# #-------# #-------# 

plot(data$behav[[1]][[index]]$position[[1]][,1], data$behav[[1]][[index]]$position[[1]][,2], 
     xlab = "x coord", ylab = "y coord",
     type="l")

# Position mouse 
position = data$behav[[1]][[index]]$position[[1]]
time_position = data$behav[[1]][[index]]$time[[1]]

saveRDS(position, file="../Data/M3424F_coord_mouse_position.RDS")
saveRDS(time_position, file="../Data/M3424F_time_mouse_position.RDS")



#-------# #-------# #-------# #-------# #-------# #-------# #-------# 
#-------# #-------#          Calcium data         #-------# #-------# 
#-------# #-------# #-------# #-------# #-------# #-------# #-------# 

neuron_data = data$neuron_dat[[1]][[index]]
rm(data)

plot(neuron_data$centroid[[1]][,1], neuron_data$centroid[[1]][,2])

loc_neurons = neuron_data$centroid[[1]]

saveRDS(loc_neurons, file="../Data/M3424F_loc_neurons.RDS")



plot(neuron_data$C_raw[[1]][1,]/max(neuron_data$C_raw[[1]][1,]), type="l",
     ylim = c(0,20))
for(i in 2:20)
{
  lines(1:length(neuron_data$C_raw[[1]][i,]),
        neuron_data$C_raw[[1]][i,]/max(neuron_data$C_raw[[1]][i,]) + i-1, 
        type="l",
       ylim = c(0,20))
}

calcium = neuron_data$C_raw[[1]]
time_calcium = neuron_data$time[[1]]

saveRDS(calcium, file="../Data/M3424F_calcium.RDS")
saveRDS(time_calcium, file="../Data/M3424F_time_calcium.RDS")

