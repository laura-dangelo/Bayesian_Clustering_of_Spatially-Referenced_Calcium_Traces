#---------# #-----------# #---------# #---------# #-----------# #---------# #---------# 
#---------#       PLOT THE DATA: CALCIUM TRACES AND MOUSE'S MOVEMENTS       #---------#
#---------# #-----------# #---------# #---------# #-----------# #---------# #---------# 
## This script is used to plot the mouse's movements and example calcium traces reported in Figure 1 of the paper

library(ggplot2)
library(ggforce)
library(ggpubr)
library(reshape2)

loc_neurons = readRDS("../Data/M3424F_loc_neurons.RDS")
data = readRDS("../Data/data_binary_position.RDS")
str(data)

idx = readRDS("../Data/Time_windows/indices.RDS")

data$pos_binary[data$pos_binary>2] = 1
which(diff(data$pos_binary)!=0)
win = unlist(apply(
  cbind(
    1:length(which(diff(data$pos_binary)!=0)),
    c(0,(which(diff(data$pos_binary)!=0)))[1:137],
    (which(diff(data$pos_binary)!=0))), 1, function(x) rep(x[1], x[3]-x[2])))

win = c(win, rep(max(win)+1, nrow(data) - length(win)))
data$win = win
str(data)



#--------------------# #--------------------# 
##        Enviroment and movements         ##
#--------------------# #--------------------# 

idcirc = (data$pos1^2 + data$pos2^2)<1.05

CIRC = ggplot()+
  theme_bw()+
  geom_circle(aes(x0=0, y0=0, r=1.03), col=1, alpha=.1,fill="royalblue1")+
  geom_circle(aes(x0=0, y0=0, r=.71), col="white", alpha=1, fill="white")+
  geom_circle(aes(x0=0, y0=0, r=.71), col=NA, alpha=.1, fill="turquoise")+
  geom_path(aes(x = data[idcirc,]$pos1, y= data[idcirc,]$pos2), alpha=.25, lwd=0.3)+
  geom_path(aes(x = data[win==17,]$pos1, y= data[win==17,]$pos2),alpha=1, col = "coral3", lwd=0.7)+
  ylab("y coordinate")+xlab("x coordinate")+
  theme_minimal() +
  theme(
    aspect.ratio=1,
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
  theme(text=element_text(size=17.5)) 
CIRC

ggsave("02_data_analysis/01_bSCDC_individual_trials/output_images/mouse_positionppp.pdf", CIRC, width = 5, height = 5)
ggsave("02_data_analysis/01_bSCDC_individual_trials/output_images/mouse_position.epf", CIRC, width = 5, height = 5, device = cairo_ps)





#--------------------# #--------------------# 
##             Calcium traces              ##
#--------------------# #--------------------# 

calcium = data[,c(16:35)]
calcium = apply(calcium, 2,function(x) x/150) #### (x-min(x))/(2+max(x)-min(x)))
TT = nrow(calcium)
n = ncol(calcium)
calcium = calcium + matrix(rep(1:n, TT), TT, n, byrow=T)
maxcl = max(calcium)

calcium = reshape2::melt((calcium))
calcium$time = rep(data$time_calcium/1000, n)
calcium$neuron = as.factor(as.numeric(calcium$Var2))

stacchi = c(1,which(abs(diff(data$pos_binary))>0),TT)
background = data.frame(
  "start" = data$time_calcium[stacchi[1:(length(stacchi)-1)]]/1000,
  "end" = data$time_calcium[stacchi[2:(length(stacchi))]]/1000,
  "y_min" = rep(0, length(stacchi)-1), 
  "y_max" = rep(maxcl, length(stacchi)-1), 
  "coll" = rep(1:2, 70)[1:(length(stacchi)-1)]
)

b17 = background[17,]  
background = background[c(1:16,18:138),]

p1 = ggplot() + 
  geom_line(data = calcium, aes(x = time, y = value, group = neuron))+
  geom_rect(data = background, xmin = background$start, xmax = background$end,
            ymin = background$y_min, ymax = background$y_max,
            fill = c("royalblue1", "lightblue")[background$coll], col = NA, alpha = 0.15)+
  geom_rect(data = b17, xmin = b17$start, xmax = b17$end,
            ymin = b17$y_min, ymax = b17$y_max, fill="coral", alpha=0.5)+
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
  theme(text=element_text(size=17)) +
  xlab("Time (seconds)")  +
  ylab("Neuron") 
p1


G = ggarrange(CIRC, p1, widths = c(1,2),
          ncol=2, nrow=1,
          common.legend = TRUE, legend="bottom")
G


ggsave("02_data_analysis/01_bSCDC_individual_trials/output_images/mouse_position_and_calcium.pdf", G, width = 12, height = 5)
ggsave("02_data_analysis/01_bSCDC_individual_trials/output_images/mouse_position_and_calcium.epf", G, width = 12, height = 5, device = cairo_ps)
