
#-------------------# #-------------------# #-------------------#
##                  PLOT QUANTITIES OF INTEREST                ##
#-------------------# #-------------------# #-------------------#

# This script replicates the sensitivity study reported in Section S4.2 of the Supplementary Material.
# Specifically, this script reproduces Figure S22 and S23 of the Supplementary Material.

library(viridisLite)
library(ggplot2)
library(ggpubr)

replications_sim = 50
errors_spike_detection= readRDS( file = "03_simulation_study/01_sensitivity_study/output_RDS/errors_spike_detection.RDS")
false_positive_rates= readRDS( file = "03_simulation_study/01_sensitivity_study/output_RDS/false_positive_rates.RDS")
false_negative_rates= readRDS( file = "03_simulation_study/01_sensitivity_study/output_RDS/false_negative_rates.RDS")
aris_neuron= readRDS( file = "03_simulation_study/01_sensitivity_study/output_RDS/ARIs.RDS")


alow_seq = factor(rep(c(0,0.5,1), replications_sim))
params_seq = factor(rep(c("a301", "b41", "c101"), each = replications_sim))

dferror = data.frame("error" = c(errors_spike_detection), "a_low" = alow_seq, "gamma_params" = params_seq)
dffn = data.frame("error" = c(false_negative_rates), "a_low" = alow_seq, "gamma_params" = params_seq)
dffp = data.frame("error" = c(false_positive_rates), "a_low" = alow_seq, "gamma_params" = params_seq)

errors_full = data.frame(rbind(dffp, dffn, dferror))
colnames(errors_full)[1] = "value"
errors_full$Error = c(cbind(rep("False positive", nrow(dffp)),rep("False negative", nrow(dffp)),rep("Misclassification", nrow(dffp)) ))


p1 = ggplot(errors_full, aes(y = value, x = Error)) +
  geom_boxplot(aes(fill=a_low),  alpha = 0.3) +
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
    # legend.title=element_blank(),
    legend.text = element_text(size=12),
    strip.text = element_text(size=12),
    text = element_text(size = 15),
    axis.title.y=element_text(size=12),
    strip.background = element_rect( fill=NA, color="gray" )
  )+
  scale_fill_manual(expression(paste(bar(a), "  ")), values = c(magma(6)[c(2,4,6)]) ) +
  facet_wrap(gamma_params ~ ., ncol= 1,
             labeller = labeller(gamma_params = c("a301" = "Gamma(3, 0.1)",
                                                  "b41" = "Gamma(4, 1)",
                                                  "c101" = "Gamma(10, 1)") ) )+
  xlab("Error type")  +
  ylab("Error rate")
p1




#----------#  ARI
alow_seq = factor(rep(c(0,0.5,1), replications_sim))
params_seq = factor(rep(c("a301", "b41", "c101"), each = replications_sim))

dfari = data.frame("ARI" = c(aris_neuron), "a_low" = alow_seq, "gamma_params" = params_seq)

p2 = ggplot(dfari, aes(y = ARI, x = a_low)) +
  geom_boxplot(aes(fill=a_low), alpha = 0.3) +
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
    # legend.title=element_blank(),
    legend.text = element_text(size=12),
    strip.text = element_text(size=12),
    text = element_text(size = 15),
    axis.title.y=element_text(size=12),
    strip.background = element_rect( fill=NA, color="gray" )
  )+
  scale_fill_manual(expression(paste(bar(a), "  ")), values = c(magma(6)[c(2,4,6)]) ) +
  facet_wrap(gamma_params ~ ., ncol= 1,
             labeller = labeller(gamma_params = c("a301" = "Gamma(3, 0.1)",
                                                  "b41" = "Gamma(4, 1)",
                                                  "c101" = "Gamma(10, 1)") ) )+
  xlab(expression(bar(a)))  +
  ylab("Adjusted Rand index")
p2


G = ggarrange(p1,p2, widths = c(2,1),
              ncol=2, nrow=1,
              common.legend = TRUE, legend="bottom")
G
ggsave("03_simulation_study/01_sensitivity_study/output_images/simulation_alow_results.pdf", G, width = 8, height = 8)
ggsave("03_simulation_study/01_sensitivity_study/output_images/simulation_alow_results..eps",G, width = 8, height = 8, device = cairo_ps)





example_data = readRDS("03_simulation_study/01_sensitivity_study/simulated_data/simulated_data1.RDS")

df = data.frame("loc1" = example_data$loc[,1], "loc2" = example_data$loc[,2], 
                "Cluster" = factor(example_data$cluster_neurons))

G2 = ggplot(df, aes(x = loc1, y = loc2, fill = Cluster)) +
  geom_point(size = 3, alpha = 0.9, shape = 21,  colour = "black",) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  scale_fill_manual(
    values = c(
      "1" = "#66C2A5",  
      "2" = "#FC8D62",  
      "3" = "#8DA0CB"   
    ),
    name = "Cluster"
  ) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill='transparent', color=NA), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.border = element_rect(color = "darkgray ", fill=NA),
    axis.line.x.bottom = element_line(color="gray"),
    legend.position = "bottom",
    legend.background = element_rect(fill='transparent', color = 'transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent', color = 'transparent'), #transparent legend panel
    legend.text = element_text(size=10),
    strip.text = element_text(size=10),
    strip.background = element_rect( fill=NA, color="gray76" ),
    text = element_text(size=12)
  )+
  xlab("")  +
  ylab("") 

G2
ggsave("03_simulation_study/01_sensitivity_study/output_images/simulation_example_data.pdf", G2, width = 5, height = 5)
ggsave("03_simulation_study/01_sensitivity_study/output_images/simulation_example_data.epf", G2, width = 5, height = 5, device = cairo_ps)
