#-------# #-------# #-------# #-------# #-------# #-------# #-------# #-------# 
#                                                                             #
#                  Import original MAT file and create RDS                    #
#                                                                             #
#-------# #-------# #-------# #-------# #-------# #-------# #-------# #-------# 

# This script imports the original .mat data and saves to the Data folder the file M3424F_trial4.rds

#  _________________________________
#  YOU CAN AVOID RUNNING THIS SCRIPT 
#  ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
# In the Google Drive folder (https://drive.google.com/drive/folders/1-1xf57mZBc1usA-iCZGkp4KPF8oX_5mV?usp=sharing)
# the RDS files produced by this script are available in the "Data" folder.
# You can download these files and copy them in the Data folder of this repository.

# FIRST - RUN THE ASSOCIATED PYTHON CODE

library(rhdf5)
mice <-  c("M3411", "M3412",  "M3421F", "M3422F", "M3424F", "M3425F")   # fill in your mouse folder names

# Bypasses h5read()'s H5Lexists() pre-check, which errors on
# the case-colliding sibling names MATLAB's HDF5 writer creates
read_dataset <- function(file, path) {
  h5f <- H5Fopen(file)
  on.exit(H5Fclose(h5f), add = TRUE)
  parts <- strsplit(sub("^/", "", path), "/")[[1]]
  obj <- h5f
  opened <- list()
  for (p in parts) {
    obj <- H5Oopen(obj, p)
    opened[[length(opened) + 1]] <- obj
  }
  data <- H5Dread(obj)
  for (o in rev(opened)) H5Oclose(o)
  data
}

lookup <- read.csv("01_data_preprocessing/trial_order_lookup.csv", stringsAsFactors = FALSE)

all_trials <- list()

for (top in 1:length(mice)) {
  
  mouse <- mice[top]
  ############################## MAKE SURE TO EXTRACT AND PLACE THIS FOLDER FROM THE DATA REPO HERE:
  file <- paste0("Fig_1_3_tri_cir_sqr/", mouse)
  behav_file  <- paste0(file, "/behav.mat")
  neuron_file <- paste0(file, "/neuronIndividuals_new.mat")
  
  mouse_lookup <- lookup[lookup$mouse == mouse, ]
  mouse_lookup <- mouse_lookup[order(mouse_lookup$trial), ]
  
  trials <- lapply(1:6, function(i) {
    bg <- mouse_lookup$behav_group[i]
    ng <- mouse_lookup$neuron_group[i]
    list(
      behav = list(
        position = read_dataset(behav_file,  paste0("/#refs#/", bg, "/position")),
        time     = read_dataset(behav_file,  paste0("/#refs#/", bg, "/time"))
      ),
      neuron = list(
        A        = read_dataset(neuron_file, paste0("/#refs#/", ng, "/A")),
        C        = read_dataset(neuron_file, paste0("/#refs#/", ng, "/C")),
        C_raw    = read_dataset(neuron_file, paste0("/#refs#/", ng, "/C_raw")),
        S        = read_dataset(neuron_file, paste0("/#refs#/", ng, "/S")),
        Cn       = read_dataset(neuron_file, paste0("/#refs#/", ng, "/Cn")),
        centroid = read_dataset(neuron_file, paste0("/#refs#/", ng, "/centroid")),
        time     = read_dataset(neuron_file, paste0("/#refs#/", ng, "/time"))
      )
    )
  })
  names(trials) <- paste0("trial", 1:6)
  
  all_trials[[mouse]] <- trials
}


out_dir <- "01_data_preprocessing/all_trials_saved"
dir.create(file.path(out_dir, "by_trial"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "by_mouse"), recursive = TRUE, showWarnings = FALSE)

all_trials <- readRDS("mice_trials.RDS")

for (mouse in names(all_trials)) {
  
  mouse_trials <- all_trials[[mouse]]
  
  # save each trial on its own
  for (trial_name in names(mouse_trials)) {
    saveRDS(
      mouse_trials[[trial_name]],
      file.path(out_dir, "by_trial", paste0(mouse, "_", trial_name, ".rds"))
    )
    message(mouse, " ", trial_name, " saved")
  }
  
  # save all 6 trials for this mouse together
  saveRDS(
    mouse_trials,
    file.path(out_dir, "by_mouse", paste0(mouse, ".rds"))
  )
  
  message("=== ", mouse, " complete ===")
}


