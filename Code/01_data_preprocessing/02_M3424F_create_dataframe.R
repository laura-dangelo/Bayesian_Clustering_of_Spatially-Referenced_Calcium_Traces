
# This script standardizes the coordinates of the mouse position and creates the two areas (inner circle and outer ring).
# It imports the individual quantities and organize them into a single, coherent, dataframe.
# It then creates the individual "trials" defined by the mouse position.
# These trials are saved into the folder "Data/Time_windows".


position = readRDS(file="../Data/M3424F_coord_mouse_position.RDS")
time_position = readRDS(file="../Data/M3424F_time_mouse_position.RDS")
calcium = readRDS(file="../Data/M3424F_calcium.RDS")
time_calcium = readRDS(file="../Data/M3424F_time_calcium.RDS")
loc_neurons = readRDS(file="../Data/M3424F_loc_neurons.RDS")


calcium_df = data.frame("time" = time_calcium[seq(1, length(time_calcium)-1, by=2)],
                        t(calcium))
rm(calcium)
rm(time_calcium)


position_df = data.frame("time" = time_position, "position1" = position[,1], "position2" = position[,2])
rm(position)
rm(time_position)

position_df = position_df[seq(1,nrow(position_df),by=2),]
position_df = position_df[1:nrow(calcium_df),]

data = data.frame("pos1" = position_df$position1,
                  "pos2" = position_df$position2,
                  "time_pos" = position_df$time,
                  "time_calcium" = calcium_df$time,
                  calcium_df[,-c(1)]
                  )
str(data)

plot(data$time_calcium, data[,5]/(max(data[,5])+max(data[,5])/4), type="l", ylim = c(0,20))
for(i in 2:20) {
  lines(data$time_calcium, data[,5+i]/(max(data[,5+i])+max(data[,5+i])/4)+i-1)
}




##-------------## ##-------------## ##-------------## ##-------------## ##-------------## ##-------------## 
## Create division between movements in the inner circle and boundary of the arena

# normalize spatial coordinates
aa = (data$pos1 - (max(data$pos1)))
bb = (data$pos2 - max(data$pos2))
aa = aa - min(aa)/2
bb = bb - min(bb)/2
aa = aa / max(aa)
bb = bb / max(bb)
plot(aa, bb, type="l")
data$pos1 = aa
data$pos2 = bb


# Normalize the coordinates so that the center of the environment is (0,0)
# and the radius is 1.
# We divide the center and the border so that the area of the two sections is similar.
# Total area is 3.14: if the radius is 0.7069 the two parts have approximately equal area.

r = 0.7069

tmp = data$pos1^2 + data$pos2^2
pos = rep(2, length(data$pos1))
pos[tmp < r^2] = 1
table(pos) # balanced observations in the two areas


data = data.frame(pos, data)
colnames(data)[1] = "pos_binary"
str(data)


# # split window 107 which is over twice as long as the others

data$pos_binary[7867:(7867+256)] = 3
data$pos_binary[(7867+256):8378] = 4

tmp = which(diff(data$pos_binary)!=0)
tmp = c(0,tmp)
tmp = c(tmp, length(data$pos1))

idx = which(diff(tmp)>50)

saveRDS(idx, file = "../Data/Time_windows/indices.RDS")

saveRDS(data, file="../Data/data_binary_position.RDS") # full data frame



##-------------## ##-------------## ##-------------## ##-------------## ##-------------## ##-------------## 
# save individual windows

for(n_window in idx )
{
  calcium = data[(tmp[n_window]+1):tmp[n_window+1] ,c(6:(ncol(data)))]
  filename = paste0("../Data/Time_windows/calcium_window", n_window, ".RDS")
  saveRDS(calcium, file = filename)
}

for(n_window in idx )
{
  time_calcium = data[(tmp[n_window]+1):tmp[n_window+1] ,5]
  filename = paste0("../Data/Time_windows/time_calcium_window", n_window, ".RDS")
  saveRDS(time_calcium, file = filename)
}

