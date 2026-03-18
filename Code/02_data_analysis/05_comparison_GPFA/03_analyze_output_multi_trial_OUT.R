library(ggplot2)
library(ggthemes)

filename = paste0("02_data_analysis/05_comparison_GPFA/RUN_GPFA_multi_trial_OUT/scale_factors_OUT.csv")
scale_fac = read.csv(file = filename)
str(scale_fac)

# study relevance of factors
scale_fac = data.frame("logs" = scale_fac$dim_scales, "nu" = (scale_fac$nus))
scale_fac$color = rank(apply( cbind(scale(scale_fac[,1]), scale(scale_fac[,2]))  ,1,mean))

p1 = ggplot(scale_fac, aes(x = logs, y = nu, fill = color)) +
  geom_point(size=3, colour="black",pch=21, alpha = 0.7) +
  xlab(expression("log s"["d"])) + 
  ylab("Latent mean scale") +
  scale_fill_gradient2_tableau(
    palette = "Sunset-Sunrise Diverging" ) +
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
    legend.direction='vertical',
    legend.position = "none",
    legend.background = element_rect(fill='transparent', color = 'white'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent'), #transparent legend panel
    legend.text = element_text(size=13),
    strip.text = element_text(size=13),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    strip.background = element_rect( fill=NA, color="gray" )
  )
p1
ggsave(p1, file = "02_data_analysis/05_comparison_GPFA/output_images/scale_vs_nu_OUT.pdf",
       width = 5, height = 4)

sort(scale_fac$logs, index.return=T)$ix
sort(scale_fac$nu, index.return=T)$ix

which((scale_fac$logs > -4.195)&(scale_fac$nu > 0.000015))
which(scale_fac$color>8)

filename = paste0("02_data_analysis/05_comparison_GPFA/RUN_GPFA_multi_trial_OUT/output_GPFA_OUT.csv")
out = read.csv(file = filename)
str(out)


df_GP = data.frame("y" = c(out[,17]), "fac" = rep(c("1"),each=45), time = rep(1:45, 2)/15)

plot_GPs = ggplot(data = df_GP, aes(x = time, y = y)) +
  geom_line(aes(y = y), col = "purple", lwd = 1 ) +
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
    legend.position = "bottom",
    legend.background = element_rect(fill='transparent', color = 'white'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent'), #transparent legend panel
    legend.text = element_text(size=10),
    strip.text = element_text(size=10),
    strip.background = element_rect( fill=NA, color="gray" )
  )+
  # scale_y_continuous(breaks = c(0.0, 0.1, 0.2,0.3),  limits = c(0,0.3)) +
  ylab("Latent factor")  + xlab("Time (seconds)") +
  facet_wrap((fac)~., ncol=1, scales = "free")
plot_GPs

