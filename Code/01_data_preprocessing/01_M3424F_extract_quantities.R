#-------# #-------# #-------# #-------# #-------# #-------# 
#                                                         #
#        Import M3424F_trial4.RDS                         #
#        and extract individual quantities                #
#                                                         #
#-------# #-------# #-------# #-------# #-------# #-------# 

# This script imports the file M3424F_trial4.RDS, extracts the individual quantities of interest and
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



M3424F_trial4 = readRDS("../Data/M3424F_trial4.rds")

str(M3424F_trial4)

loc_neurons = M3424F_trial4$neuron$centroid
plot(loc_neurons[,1], loc_neurons[,2])
saveRDS(loc_neurons, file = "../Data/M3424F_loc_neurons.RDS")

time_calcium = M3424F_trial4$neuron$time
saveRDS(time_calcium, file = "../Data/M3424F_time_calcium.RDS")

time_mouse_position = M3424F_trial4$behav$time
saveRDS(time_mouse_position, file = "../Data/M3424F_time_mouse_position.RDS")

calcium = M3424F_trial4$neuron$C_raw
saveRDS(calcium, file = "../Data/M3424F_calcium.RDS")

coord_mouse = M3424F_trial4$behav$position
saveRDS(coord_mouse, file = "../Data/M3424F_coord_mouse_position.RDS")
