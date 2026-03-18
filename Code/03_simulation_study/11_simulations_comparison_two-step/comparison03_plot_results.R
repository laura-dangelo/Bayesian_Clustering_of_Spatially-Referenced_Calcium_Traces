
errors_spike_detection_JW= readRDS( file = "11_simulations_comparison_two-step/output_RDS/errors_spike_detection_JW.RDS")
false_positive_rates_JW= readRDS( file = "11_simulations_comparison_two-step/output_RDS/false_positive_rates_JW.RDS")
false_negative_rates_JW= readRDS( file = "11_simulations_comparison_two-step/output_RDS/false_negative_rates_JW.RDS")


errors_spike_detection= readRDS( file = "11_simulations_comparison_two-step/output_RDS/errors_spike_detection.RDS")
false_positive_rates= readRDS( file = "11_simulations_comparison_two-step/output_RDS/false_positive_rates.RDS")
false_negative_rates= readRDS( file = "11_simulations_comparison_two-step/output_RDS/false_negative_rates.RDS")

errors_spike_detection = errors_spike_detection[,9]
false_positive_rates = false_positive_rates[,9]
false_negative_rates = false_negative_rates[,9]


#----------# #----------# #----------# #----------#
#----------# PLOT ERROR RATE OF SPIKE DETECTION

library(ggplot2)
Fpar = c(rep("Simultaneous deconvolution \n & spatial clustering", 50), rep("L0 + Consensus K-Means", 50))
dferror = data.frame("error" = c(errors_spike_detection, errors_spike_detection_JW), "par" = Fpar)
dferror$par = as.factor(dferror$par)


library(viridisLite)
errors_full = data.frame("value" = c(c(false_positive_rates), c(false_negative_rates),c(errors_spike_detection),
                                     c(c(false_positive_rates_JW), c(false_negative_rates_JW),c(errors_spike_detection_JW))),
                         "Error" = c(rep("FP", 50),rep("FN", 50),rep("Misclass.", 50) ,
                                     rep("FP", 50),rep("FN", 50),rep("Misclass.", 50) ),
                         "Model" = c(rep("Simultaneous deconvolution \n & spatial clustering", 150), rep("L0 + Consensus KMeans", 150))
)


p1 = ggplot(errors_full, aes(y = value, x = Error)) +
  geom_boxplot(aes( fill=Model),  alpha = 0.3) +
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
    legend.background = element_rect(fill='transparent', color = 'transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent', color = 'transparent'), #transparent legend panel
    legend.title=element_blank(),
    legend.text = element_text(size=12),
    strip.text = element_text(size=12),
    text = element_text(size = 15),
    strip.background = element_rect( fill=NA, color="gray" )
  )+
  # scale_fill_manual( values = c(rocket(8)[4], mako(8)[5], inferno(8)[7]) ) +
  # scale_fill_manual( values = c("cyan4", "deeppink4", "darkgoldenrod1") ) +
  scale_fill_manual( values = c(magma(6)[c(2,4)]) ) +
  ylim(c(0,1)) +
  xlab("")  +
  ylab("Error rate")
p1


#----------# #----------# #----------# #----------#
#----------# PLOT ARI

name = paste0("11_simulations_comparison_two-step/output_RDS/ARIs_kmeans.RDS")
ARIs_kmeans = readRDS(file = name)

name = paste0("11_simulations_comparison_two-step/output_RDS/ARIs.RDS")
ARIs = readRDS(file = name)
ARIs = ARIs[,9]


#----------#  ARI
dferror = data.frame("error" = c(ARIs, ARIs_kmeans), "Model" = Fpar)
dferror$Model = as.factor(dferror$Model)

p2 = ggplot(dferror, aes(y = error, x = Model, fill=Model)) +
  geom_boxplot(  alpha = 0.3) +
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
    legend.background = element_rect(fill='transparent', color = 'transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent', color = 'transparent'), #transparent legend panel
    legend.title=element_blank(),
    legend.text = element_text(size=13),
    strip.text = element_text(size=12),
    text = element_text(size = 15),
    axis.title.y=element_text(size=12),
    strip.background = element_rect( fill=NA, color="gray" ),
    # axis.text.x = element_blank(),
    # axis.text.y = element_blank(),
    axis.ticks = element_blank()
  )+
  # scale_fill_manual( values = c(rocket(8)[4], mako(8)[5], inferno(8)[7]) ) +
  # scale_fill_manual( values = c("cyan4", "deeppink4", "darkgoldenrod1") ) +
  scale_fill_manual( values = c(magma(6)[c(2,4)]) ) +
  scale_x_discrete(labels=c("","")) +
  ylim(c(0,1)) +
  xlab("")  +
  ylab("ARI")
p2



library(ggpubr)
G = ggarrange(p1,p2, widths = c(2,1),
              ncol=2, nrow=1,
              common.legend = TRUE, legend="bottom")
G
ggsave("11_simulations_comparison_two-step/output_images/simulation_results.pdf", G, width = 8, height = 4)
ggsave("11_simulations_comparison_two-step/output_images/simulation_results.eps",G, width = 8, height = 4, device = cairo_ps)

