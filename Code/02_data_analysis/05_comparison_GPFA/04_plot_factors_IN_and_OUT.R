
library(patchwork)
library(ggplot2)

filename = paste0("02_data_analysis/05_comparison_GPFA/RUN_GPFA_multi_trial_OUT/output_GPFA_OUT.csv")
out_OUT = read.csv(file = filename)
str(out_OUT)

filename = paste0("02_data_analysis/05_comparison_GPFA/RUN_GPFA_multi_trial_IN/output_GPFA_IN.csv")
out_IN = read.csv(file = filename)
str(out_IN)


df_GP = data.frame("y" = c(out_IN[,1],out_IN[,2],out_IN[,24], out_OUT[,17]), 
                   "fac" = rep(c("Factor 1 - IN","Factor 2 - IN","Factor 3 - IN","Factor 1 - OUT"),each=45), time = rep(1:45, 4)/15)

df_GP$y = df_GP$y * 1000


p_GP1 = ggplot(data = df_GP[df_GP$fac == "Factor 1 - IN",], aes(x = time, y = y)) +
  geom_line(aes(y = y), col = "royalblue3", lwd = 1 ) +
  theme_minimal() +
  theme(
    panel.spacing = unit(1, "lines"),
    panel.background = element_rect(fill='transparent', color=NA), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    # panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    panel.border = element_rect(color = "darkgray ", fill=NA),
    # axis.line.y.left = element_line(color="gray"),
    axis.line.x.bottom = element_line(color="gray"),
    #
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom",
    legend.background = element_rect(fill='transparent', color = 'white'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent'), #transparent legend panel
    legend.text = element_text(size=10),
    strip.text = element_text(size=10),
    strip.background = element_rect( fill=NA, color="gray" )
  )+
  ylab("Latent factor")  + xlab("Time (seconds)") + ggtitle("Factor 1 - IN")
p_GP1

p_GP2 = ggplot(data = df_GP[df_GP$fac == "Factor 2 - IN",], aes(x = time, y = y)) +
  geom_line(aes(y = y), col = "royalblue3", lwd = 1 ) +
  theme_minimal() +
  theme(
    panel.spacing = unit(1, "lines"),
    panel.background = element_rect(fill='transparent', color=NA), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    # panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    panel.border = element_rect(color = "darkgray ", fill=NA),
    # axis.line.y.left = element_line(color="gray"),
    axis.line.x.bottom = element_line(color="gray"),
    #
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom",
    legend.background = element_rect(fill='transparent', color = 'white'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent'), #transparent legend panel
    legend.text = element_text(size=10),
    strip.text = element_text(size=10),
    strip.background = element_rect( fill=NA, color="gray" )
  )+
  ylab(" ")  + xlab("Time (seconds)") + ggtitle("Factor 2 - IN")
p_GP2

p_GP3 = ggplot(data = df_GP[df_GP$fac == "Factor 3 - IN",], aes(x = time, y = y)) +
  geom_line(aes(y = y), col = "royalblue3", lwd = 1 ) +
  theme_minimal() +
  theme(
    panel.spacing = unit(1, "lines"),
    panel.background = element_rect(fill='transparent', color=NA), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    # panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    panel.border = element_rect(color = "darkgray ", fill=NA),
    # axis.line.y.left = element_line(color="gray"),
    axis.line.x.bottom = element_line(color="gray"),
    #
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom",
    legend.background = element_rect(fill='transparent', color = 'white'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent'), #transparent legend panel
    legend.text = element_text(size=10),
    strip.text = element_text(size=10),
    strip.background = element_rect( fill=NA, color="gray" )
  )+
  ylab(" ")  + xlab("Time (seconds)") + ggtitle("Factor 3 - IN")
p_GP3


p_GP4 = ggplot(data = df_GP[df_GP$fac == "Factor 1 - OUT",], aes(x = time, y = y)) +
  geom_line(aes(y = y), col = "orchid4", lwd = 1 ) +
  theme_minimal() +
  theme(
    panel.spacing = unit(1, "lines"),
    panel.background = element_rect(fill='transparent', color=NA), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    # panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    panel.border = element_rect(color = "darkgray ", fill=NA),
    # axis.line.y.left = element_line(color="gray"),
    axis.line.x.bottom = element_line(color="gray"),
    #
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom",
    legend.background = element_rect(fill='transparent', color = 'white'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent'), #transparent legend panel
    legend.text = element_text(size=10),
    strip.text = element_text(size=10),
    strip.background = element_rect( fill=NA, color="gray" )
  )+
  ylab("Latent factor")  + xlab("Time (seconds)") + ggtitle("Factor 1 - OUT")
p_GP4



pp=wrap_plots(
  p_GP1, p_GP2, p_GP3, 
  ncol = 3
)
pp
ggsave(pp, file = "02_data_analysis/05_comparison_GPFA/output_images/lGP_time_IN.pdf",
       width = 11, height =3)

p_GP4
ggsave(p_GP4, file = "02_data_analysis/05_comparison_GPFA/output_images/lGP_time_OUT.pdf",
       width = 11/3, height =3)
