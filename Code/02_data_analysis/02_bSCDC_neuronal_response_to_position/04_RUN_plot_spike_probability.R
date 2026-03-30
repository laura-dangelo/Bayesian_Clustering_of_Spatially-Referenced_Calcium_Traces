
#---------# #-----------# #---------# #---------# #-----------# #---------# 
#---------#        PLOT THE NEURONS' SPIKE PROBABILITY          #---------#
#---------# #-----------# #---------# #---------# #-----------# #---------# 
# This script produces the plots of Figure S3 of the Supplementary Material, which show how the spike probability 
# of individual neurons is affected by the mouse's position in the arena.

library(tidyverse)

data = readRDS("../Data/data_binary_position.RDS")
idx = readRDS("../Data/Time_windows/indices.RDS")

data$pos_binary[7867:(7867+256)] = 3
data$pos_binary[(7867+256):8378] = 4

tmp = which(diff(data$pos_binary)!=0)
tmp = c(0,tmp)
tmp = c(tmp, length(data$pos1))

idx = which(diff(tmp)>50)

in_or_out <- matrix(NA,length(idx),2)
colnames(in_or_out) <- c("index","in_or_out")
o = 0
for(n_window in idx )
{
  o = o+1
  in_or_out[o,1] <- n_window
  in_or_out[o,2] <- unique(data[(tmp[n_window]+1):tmp[n_window+1],1])
}
in_or_out = in_or_out %>% as_tibble()


locat <- readRDS("../Data/M3424F_loc_neurons.RDS") 


in_  <- in_or_out %>% filter(in_or_out == 1)
out_ <- in_or_out %>% filter(in_or_out == 2)

SPin <- CL_in <- array(NA,c(325,nrow(in_)))

for(i in 1:nrow(in_)){
  
  SPin[,i]  <- colMeans(readRDS(paste0("02_data_analysis/01_bSCDC_individual_trials/output_RDS/estimated_spikes_win",
                                       in_[i,]$index,".RDS")))

}


SPout <- CL_out <- array(NA,c(325,nrow(out_)))

for(i in 1:nrow(out_)){
  
  SPout[,i]  <- colMeans(readRDS(paste0("02_data_analysis/01_bSCDC_individual_trials/output_RDS/estimated_spikes_win",
                                        out_[i,]$index,".RDS")))
}




layout_matrix_in <- data.frame(locat)
layout_matrix_in = cbind(layout_matrix_in, "var" = rep("Inner circle", 325), "col" = rowMeans(SPin))
layout_matrix_in = layout_matrix_in[order(layout_matrix_in$col),]

layout_matrix_out <- data.frame(locat)
layout_matrix_out = cbind(layout_matrix_out, "var" = rep("Outer ring", 325), "col" = rowMeans(SPout))
layout_matrix_out = layout_matrix_out[order(layout_matrix_out$col),]

layout_matrix = rbind(layout_matrix_in, layout_matrix_out)


pnum = ggplot() + 
  geom_point( aes(x = X1, y = X2, 
                  fill  = col,
                  size  = col*0.8), 
              data = layout_matrix, 
              alpha=.9, col=1, pch=21) +
  scale_fill_gradient2(
    low = "black",
    mid = "deeppink2",
    high = "yellow",
    midpoint = 0.08,
    # limits = c(0,0.15), 
    name = "Spike probability"
  ) +
  guides(size = "none") +
  theme_bw() +
  theme(
    aspect.ratio=1,
    panel.background = element_rect(fill='transparent', color=NA), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.minor = element_blank(), #remove minor gridlines
    panel.border = element_rect(color = "darkgray ", fill=NA),
    axis.line.x.bottom = element_line(color="gray"),
    #
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5),
    legend.title = element_text(vjust = 0.75, size = 11),
    legend.background = element_rect(fill='transparent', color = 'transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent', color = 'transparent'), #transparent legend panel
    legend.text = element_text(size=10),
    strip.text = element_text(size=11),
    title = element_text(size=11),
    strip.background = element_rect( fill=NA, color="gray76" ),
    axis.title=element_text(size=11))+
  xlab("x coordinate") + 
  ylab("y coordinate") + 
  facet_wrap(~var) 
pnum

ggsave(pnum, file="02_data_analysis/02_bSCDC_neuronal_response_to_position/output_images/Average_spikes.pdf", height = 5.5, width = 9)



