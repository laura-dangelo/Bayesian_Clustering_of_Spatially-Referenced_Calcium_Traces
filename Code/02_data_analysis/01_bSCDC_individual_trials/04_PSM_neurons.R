
#---------# #-----------# #---------# #---------# #-----------# #---------# 
#---------#   COMPUTE AND PLOT POSTERIOR SIMILARITY MATRICES    #---------#
#---------# #-----------# #---------# #---------# #-----------# #---------# 
# This script produces the posterior similarity matrix of the neurons' functional cluster allocation,
# similar to Figure 6 in the main paper.

library(ggplot2)
library(ggforce)
library(salso)

# auxiliary function: find mode of a discrete vector
find_mode <- function(x) {
  freq_table <- table(x)
  mode <- as.numeric(names(freq_table[freq_table == max(freq_table)]))
  return(mode)
}

# import data
idx <- readRDS("../Data/Time_windows/indices.RDS")

if(!file.exists("02_data_analysis/02_bSCDC_neuronal_response_to_position/output_RDS/uncertClust.RDS")){
  PSMs = c()
  Loss = c()
  Chai = c()
  for(n_window in idx){
    
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
              chai = Chai),file = "02_data_analysis/02_bSCDC_neuronal_response_to_position/output_RDS/uncertClust.RDS")
  
}else{
 l = readRDS("02_data_analysis/02_bSCDC_neuronal_response_to_position/output_RDS/uncertClust.RDS")
 PSMs =  l$psm
 Loss =  l$loss
 Chai =  l$chai
}

rm(l)
summary(PSMs$value)
colnames(PSMs)
head(PSMs)

PSMs$name_win <- "Window"
PSMs$name_win[PSMs$window==4] <- "Window #04"
PSMs$name_win[PSMs$window==17] <- "Window #17"
PSMs$name_win[PSMs$window==36] <- "Window #36"
PSMs$name_win[PSMs$window==77] <- "Window #77"

# plot the matrices for the representative windows

psms17 = ggplot(PSMs[PSMs$window==17,])+
  # theme_void()+
  theme_minimal() +
  theme(
    panel.background = element_rect(fill='transparent', color=NA), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    # panel.border = element_rect(color = "darkgray ", fill=NA),
    # axis.line.y.left = element_line(color="gray"),
    axis.line.x.bottom = element_line(color="transparent"),
    #
    legend.position = "bottom",
    legend.background = element_rect(fill='transparent', color = 'transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent', color = 'transparent'), #transparent legend panel
    # legend.title=element_blank(),
    legend.text = element_text(size=8),
    strip.text = element_text(size=11),
    text = element_text(size = 11),
    strip.background = element_rect( fill=NA, color="transparent" )
  )+
  geom_tile(aes(x=Var1,
                y=Var2,
                fill=value))+
  xlab("Neuron") + ylab("Neuron") +
  scale_fill_viridis_c("Probability", option="inferno", begin = 0.15, end = 1, direction = -1)+
  theme(legend.position= "right")
psms17


filename = paste0("02_data_analysis/01_bSCDC_individual_trials/output_images/PSM17.pdf")
ggsave(filename, psms17, width = 5, height = 4)
filename = paste0("02_data_analysis/01_bSCDC_individual_trials/output_images/PSM17.png")
ggsave(filename, psms17, width = 5, height = 4)



psms1 = ggplot(PSMs[PSMs$window %in% c(4,36,77),])+
  # theme_void()+
  theme_minimal() +
  theme(
    panel.background = element_rect(fill='transparent', color=NA), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    # panel.border = element_rect(color = "darkgray ", fill=NA),
    # axis.line.y.left = element_line(color="gray"),
    axis.line.x.bottom = element_line(color="transparent"),
    #
    legend.position = "bottom",
    legend.background = element_rect(fill='transparent', color = 'transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent', color = 'transparent'), #transparent legend panel
    # legend.title=element_blank(),
    legend.text = element_text(size=8),
    strip.text = element_text(size=11),
    text = element_text(size = 11),
    strip.background = element_rect( fill=NA, color="transparent" )
  )+
  geom_tile(aes(x=Var1,
                y=Var2,
                fill=value))+
  xlab("Neuron") + ylab("Neuron") + 
  facet_wrap(~name_win, scales='free')+
  scale_fill_viridis_c("Probability", option="inferno", begin = 0.15, end = 1, direction = -1)+
  theme(legend.position= "bottom")
psms1



filename = paste0("02_data_analysis/01_bSCDC_individual_trials/output_images/PSM.pdf")
ggsave(filename, psms1, width = 10, height = 4.5)
filename = paste0("02_data_analysis/01_bSCDC_individual_trials/output_images/PSM.png")
ggsave(filename, psms1, width = 10, height = 4.5)



