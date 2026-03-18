
#---------# #-----------# #---------# #---------# #-----------# #---------# 
#---------#                DATA PREPARATION GPFA                #---------#
#---------# #-----------# #---------# #---------# #-----------# #---------# 
# Data preparation in the 3-d array format used for the multi-trial formulation.

library(tidyr)
library(dplyr)
library(ggplot2)
library(tidyverse)

data = readRDS("../Data/data_binary_position.RDS")
idx = readRDS("../Data/Time_windows/indices.RDS")

a = readRDS("../Data/Time_windows/time_calcium_window1.RDS")
p = readRDS("../Data/M3424F_coord_mouse_position.RDS")
tp = readRDS("../Data/M3424F_time_mouse_position.RDS")


data$pos_binary[7867:(7867+256)] = 3
data$pos_binary[(7867+256):8378] = 4

tmp = which(diff(data$pos_binary)!=0)
tmp = c(0,tmp)
tmp = c(tmp, length(data$pos1))

in_or_out <- matrix(NA,length(idx),2)
colnames(in_or_out) <- c("index","in_or_out")
o = 0
for(n_window in idx )
{
  o = o+1
  in_or_out[o,1] <- n_window
  in_or_out[o,2] <- unique(data[(tmp[n_window]+1):tmp[n_window+1],1])
}


binary_pos = data$pos_binary
n = 325
TT = nrow(data)
rm(data)
counts_spikes = matrix(0, n, TT)

for(i in 1:325){
  name = paste0("02_data_analysis/05_comparison_GPFA/results/list_spikes_JW_neuron", i, ".RDS")
  tmp = readRDS(name)
  
  name = paste0("02_data_analysis/05_comparison_GPFA/results/number_spikes_JW_neuron", i, ".RDS")
  tmp2 = readRDS(name)
  tmp2[tmp2<0]=0
  
  counts_spikes[i,tmp] = tmp2
}


idx = which(diff(tmp)>50)

for(n_window in idx )
{
  calcium = counts_spikes[, (tmp[n_window]+1):tmp[n_window+1]]
  filename = paste0("02_data_analysis/05_comparison_GPFA/results/counts_spikes_window", n_window, ".RDS")
  saveRDS(calcium, file = filename)
}


#---------# #-----------# #---------# #---------#
#           create arrays IN and OUT            #
#---------# #-----------# #---------# #---------#

inds <- (tp < max(a)) & (tp > min(a))
in1 <- in_or_out %>% as_tibble %>%  dplyr::filter(in_or_out == 1) %>% pull(index)
in2 <- in_or_out %>% as_tibble %>%  dplyr::filter(in_or_out == 2) %>% pull(index)


x <- readRDS(paste0("02_data_analysis/05_comparison_GPFA/results/counts_spikes_window", idx[1], ".RDS"))

dim(x)
# n_neurons x n_timesteps

#  n_samples x n_timesteps x n_neurons
ARR_IN_CALCIUM <- array(NA,c(length(in1), 45, nrow(x)))
ARR_OUT_CALCIUM <- array(NA,c(length(in2), 45, nrow(x)))
str(ARR_IN_CALCIUM)
str(ARR_OUT_CALCIUM)

for(i in seq_along(in1)){
  x <- readRDS(paste0("02_data_analysis/05_comparison_GPFA/results/counts_spikes_window", in1[i], ".RDS"))
  ARR_IN_CALCIUM[i,,] <- as.matrix(x[,1:45])
  
}

for(i in seq_along(in2)){
  x <- readRDS(paste0("02_data_analysis/05_comparison_GPFA/results/counts_spikes_window", in2[i], ".RDS"))
  ARR_OUT_CALCIUM[i,,] <- as.matrix(x[,1:45])
}


dim(ARR_IN_CALCIUM)
dim(ARR_OUT_CALCIUM)

saveRDS(ARR_IN_CALCIUM,"../Data/GPFA_array_in_circle_NumTrials_x_Time45_x_Neurons.RDS")
saveRDS(ARR_OUT_CALCIUM,"../Data/GPFA_array_out_circle_NumTrials_x_Time45_x_Neurons.RDS")

# in the paper we used
# 4, 17, 36, 77

for(win in c(4, 17, 36, 77)){
  if(win %in% in1) print(paste0("Window ", win, " is inside"))
  if(win %in% in2) print(paste0("Window ", win, " is outside"))
}





# CHECK
ARR_OUT_CALCIUM = readRDS(file="../Data/GPFA_array_out_circle_NumTrials_x_Time45_x_Neurons.RDS")
str(ARR_OUT_CALCIUM)

ARR_IN_CALCIUM = readRDS(file="../Data/GPFA_array_in_circle_NumTrials_x_Time45_x_Neurons.RDS")
str(ARR_IN_CALCIUM)


index_win = 3 # pick index

dat2 <-
  ARR_OUT_CALCIUM[index_win,,] %>%
  as_tibble() %>%
  rownames_to_column("Var1") %>%
  pivot_longer(-Var1, names_to = "Var2", values_to = "value") %>%
  mutate(
    Var1 = factor(Var1, levels = 1:45),
    Var2 = factor(gsub("V", "", Var2), levels = 1:325)
  )
ggplot(dat2, aes(Var1, Var2)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient(low = "white", high = "red") +
  scale_y_discrete(breaks=seq(1,325,by=10))



dat3 <-
  ARR_IN_CALCIUM[index_win,,] %>%
  as_tibble() %>%
  rownames_to_column("Var1") %>%
  pivot_longer(-Var1, names_to = "Var2", values_to = "value") %>%
  mutate(
    Var1 = factor(Var1, levels = 1:45),
    Var2 = factor(gsub("V", "", Var2), levels = 1:325)
  )
ggplot(dat3, aes(Var1, Var2)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient(low = "white", high = "red") +
  scale_y_discrete(breaks=seq(1,325,by=10))


