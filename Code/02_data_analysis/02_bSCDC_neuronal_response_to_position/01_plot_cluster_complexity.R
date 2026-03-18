
#---------# #-----------# #---------# #---------# #-----------# #---------# 
#---------#   COMPUTE AND PLOT INDICES OF CLUSTER COMPLEXITY    #---------#
#---------# #-----------# #---------# #---------# #-----------# #---------# 
# This script produces the plots of the cluster complexity in relation to the mouse's position in the arena,
# similar to Figure 8 of the main paper.
# Specifically, it shows how the variance and the mode of the number of clusters changes depending on the spatial location.


library(tidyverse)
library(dplyr)
library(tidyr)   # For expand_grid

library(ggplot2)
library(fields)  # For image.smooth
library(MBA)     # For mba.surf
library(ggforce) # For geom_circle
library(patchwork)
library(interp)
library(fields)


idx  <- readRDS("../Data/Time_windows/indices.RDS")

# If missing, create Windows_list.RDS
if(!file.exists("02_data_analysis/02_bSCDC_neuronal_response_to_position/output_RDS/Windows_list.RDS")){
  WIND <- list()
  nix <- 0
  for(n_window in idx){
    
    nix = nix+1
    filename = paste0("02_data_analysis/01_bSCDC_individual_trials/results/run_gibbs_window", n_window, "_alow_1.RDS")
    out = readRDS(file = filename)
    
    a = out$AA>0
    rm(out)
    gc()
    
    WIND[[nix]] <- t(apply(a,c(1,2),mean))
    cat(nix)
  }
  saveRDS(WIND,"02_data_analysis/02_bSCDC_neuronal_response_to_position/output_RDS/Windows_list.RDS")
} else {
  WIND <- readRDS("02_data_analysis/02_bSCDC_neuronal_response_to_position/output_RDS/Windows_list.RDS")
}


# import data
data <- readRDS("../Data/data_binary_position.RDS")

data$pos_binary[7867:(7867+256)] <- 3
data$pos_binary[(7867+256):8378] <- 4
idcirc = (data$pos1^2 + data$pos2^2)<1.05
win = unlist(apply(
  cbind(
    1:length(which(diff(data$pos_binary)!=0)),
    c(0,(which(diff(data$pos_binary)!=0)))[1:138],
    (which(diff(data$pos_binary)!=0))), 1, function(x) rep(x[1], x[3]-x[2])))
win = c(win, rep(139, 10870-10862))
data$win = win
data = data[,c(1:5,331,6:330)]


# import uncertClust.RDS (or create it, if not available)
if(!file.exists("02_data_analysis/02_bSCDC_neuronal_response_to_position/output_RDS/uncertClust.RDS")){
  PSMs = c()
  Loss = c()
  Chai = c()
  for(n_window in idx_to_run){
    
    filename = paste0("02_data_analysis/01_bSCDC_individual_trials/results/run_gibbs_window", n_window, "_alow_1.RDS")
    out_alow1 = readRDS(file = filename)
    
    a    <- (apply(t(out_alow1$cluster_signal),1,function(y) c( length(unique(y)))))
    Chai <- rbind(Chai,cbind(win = n_window, mea = mean(a), mod = find_mode(a), var = var(a)))
    psm <- salso::psm(t(out_alow1$cluster_signal),nCores = 5)
    
    # sort by estimated cluster 
    filename = paste0("02_data_analysis/01_bSCDC_individual_trials/output_RDS/est_cluster_neurons_win", n_window, ".RDS")
    est_cluster_neurons = readRDS(file=filename)
    ord = sort.int(est_cluster_neurons, index.return = TRUE)$ix
    
    x <- salso::salso(t(out_alow1$cluster_signal))
    PSMs <- rbind(PSMs, cbind(reshape2::melt(psm[ord,ord]), window = n_window))
    Loss <- rbind(Loss, cbind(win = n_window, el = attr(x,"info")[4]))
    cat(n_window)  
  }
  saveRDS(list(psm = PSMs,
               loss = Loss,
               chai = Chai), file = "02_data_analysis/02_bSCDC_neuronal_response_to_position/output_RDS/uncertClust.RDS")
  
}else{
  l = readRDS("02_data_analysis/02_bSCDC_neuronal_response_to_position/output_RDS/uncertClust.RDS")
  PSMs =  l$psm
  Loss =  l$loss
  Chai =  l$chai
}
rm(l)


unc_clust_ggplot = c()
for(i in idx){
  unc_clust_ggplot = rbind(unc_clust_ggplot,
                           cbind(x = data[win == i,]$pos1, 
                                 y = data[win == i,]$pos2,
                                 win = Chai[Loss$win == i,1],
                                 mod = Chai[Loss$win == i,3],
                                 var = Chai[Loss$win == i,4]  )
  )  
}
unc_clust_ggplot <- unc_clust_ggplot %>% as_tibble()





#--------------------# #--------------------# #--------------------# #--------------------# 
#   Plot the smoothed surface of the variance of the clustering
#--------------------# #--------------------# #--------------------# #--------------------# 

set.seed(123)
n <- nrow(unc_clust_ggplot)
x <- unc_clust_ggplot$x
y <- unc_clust_ggplot$y
z <- unc_clust_ggplot$var
data <- cbind(x, y, z)

fit <- mba.surf(
  data,
  no.X = 200,
  no.Y = 200,
  extend = FALSE
)

new <- expand_grid(fit$xyz.est$x,fit$xyz.est$y)
newx <- new[,1]
newy <- new[,2]

fit$xyz.est$z2 <- image.smooth(fit$xyz.est$z, theta = 10)$z
fit$xyz.est$z2[is.na(fit$xyz.est$z)] <- NA

V <- ggplot()+
  geom_tile(aes(x=pull(newx),
                y=pull(newy),
                fill=round(c(fit$xyz.est$z2),4)))+
  scale_fill_gradient2("Variance",low = "transparent",
                       mid = "lightpink1",
                       high = "maroon",
                       midpoint = .12,
                       na.value = "transparent",
                       breaks = seq(-0.05, 0.25, 0.05),
                       labels = c("", seq(0, 0.25, 0.05) )
                       )+
 geom_contour(color = "black", alpha = 0.1)+
  geom_circle(aes(x0=0, y0=0, r=1.05), col=1, alpha=0,fill="transparent")+
  geom_circle(aes(x0=0, y0=0, r=.73), col=1, alpha=0,fill="transparent")+
  theme_minimal() +
  guides(fill = guide_colorsteps(barwidth = unit(10, "cm"),
                                 barheight = unit(.25, "cm"),
                                 show.limits = F) )+
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
    legend.title = element_text(vjust = 0.75, size = 15),
    legend.background = element_rect(fill='transparent', color = 'transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent', color = 'transparent'), #transparent legend panel
    legend.text = element_text(size=12),
    strip.text = element_text(size=10),
    title = element_text(size=14),
    strip.background = element_rect( fill=NA, color="gray76" ),
    text = element_text(size=14),
    axis.title=element_text(size=15)
  )+
  xlab("x coordinate")  +
  ylab("y coordinate") +
  ggtitle("Posterior variance of the number of clusters") 

V
ggsave(filename = "02_data_analysis/02_bSCDC_neuronal_response_to_position/output_images/post_variance.png",
      plot = V, width = 6, height = 6)
ggsave(filename = "02_data_analysis/02_bSCDC_neuronal_response_to_position/output_images/post_variance.pdf",
       plot = V, width = 6, height = 6)






#--------------------# #--------------------# #--------------------# #--------------------# 
#   Plot the smoothed surface of the mode of the clustering
#--------------------# #--------------------# #--------------------# #--------------------# 

n <- nrow(unc_clust_ggplot)
x <- unc_clust_ggplot$x
y <- unc_clust_ggplot$y
z <- unc_clust_ggplot$mod
data <- cbind(x, y, z)

fit <- mba.surf(
  data,
  no.X = 200,
  no.Y = 200,
  extend = FALSE
)


new <- expand_grid(fit$xyz.est$x,fit$xyz.est$y)
newx <- new[,1]
newy <- new[,2]

fit$xyz.est$z2 <- image.smooth(fit$xyz.est$z, theta = 10)$z
fit$xyz.est$z2[is.na(fit$xyz.est$z)] <- NA

M <- ggplot()+
  geom_tile(aes(x=pull(newx),
                y=pull(newy),
                fill=round(c(fit$xyz.est$z2),4)))+
  scale_fill_gradient2("Mode",low = "transparent",
                       mid = "lightblue",
                       high = "darkblue",
                       midpoint = 14,
                       na.value = "transparent",
                       breaks = seq(7,20,2),
                       labels = seq(7,20,2) )+
  geom_contour(color = "black", alpha = 0.1)+
  geom_circle(aes(x0=0, y0=0, r=1.05), col=1, alpha=0,fill="transparent")+
  geom_circle(aes(x0=0, y0=0, r=.73), col=1, alpha=0,fill="transparent")+
  theme_minimal() +
  guides(fill = guide_colorsteps(barwidth = unit(10, "cm"),
                                 barheight = unit(.25, "cm"), 
                                 show.limits = F))+
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
    legend.title = element_text(vjust = 0.75, size = 15),
    legend.background = element_rect(fill='transparent', color = 'transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent', color = 'transparent'), #transparent legend panel
    legend.text = element_text(size=12),
    strip.text = element_text(size=10),
    title = element_text(size=14),
    strip.background = element_rect( fill=NA, color="gray76" ),
    text = element_text(size=14),
    axis.title=element_text(size=15)
  )+
  xlab("x coordinate")  +
  ylab("y coordinate") +
  ggtitle("Posterior mode of the number of clusters") +
  guides(colour=guide_colourbar(barwidth=8))


M
ggsave(filename = "02_data_analysis/02_bSCDC_neuronal_response_to_position/output_images/post_mode.png",
       plot = M, width = 6, height = 6)
ggsave(filename = "02_data_analysis/02_bSCDC_neuronal_response_to_position/output_images/post_mode.pdf",
       plot = M, width = 6, height = 6)






