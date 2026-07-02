
#---------# #-----------# #---------# #---------# #-----------# #---------#
#---------#        Graphical representation of the prior        #---------#
#---------# #-----------# #---------# #---------# #-----------# #---------#

# This script produces the individual graphs in Figure 2 of the main paper.


library(ggplot2)
library(mvtnorm)
source("03_simulation_study/01_sensitivity_study/00_auxiliary_functions_DONT_RUN.R")


n = 22
TT = 80


#---------# #-----------# #---------# #---------# #-----------# #---------#
#---------# #-----------#     location neurons    #-----------# #---------#
#---------# #-----------# #---------# #---------# #-----------# #---------#

set.seed(14)
locations = rbind(
  rmvnorm(7, c(4, 6), diag(2)*0.5),
  rmvnorm(7, c(3.5, 4.5), diag(2)*0.2),
  rmvnorm(8, c(6, 5), diag(2)*0.2)
)
cluster_neurons = c( rep(1, 7), rep(2, 7), rep(3, 8) )
table(cluster_neurons)

df = data.frame(x = locations[,1],
                y = locations[,2],
                cl = factor(cluster_neurons)
)
p1 = ggplot(df, aes(x=x, y=y)) + 
  geom_point(colour="black",pch=21, size=3, alpha=.5 , fill="blue") + 
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
    legend.position = "none",
    legend.background = element_rect(fill='transparent', color = 'transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent', color = 'transparent'), #transparent legend panel
    # legend.title=element_blank(),
    legend.text = element_text(size=12),
    strip.text = element_text(size=12),
    text = element_text(size = 12),
    axis.title.x=element_text(size=12),
    axis.title.y=element_text(size=12),
    strip.background = element_rect( fill=NA, color="gray" )
  ) + 
  xlab("x coordinate") + 
  ylab("y coordinate")  
# scale_fill_manual(values = c("deeppink3", "#00BFC4", "blue"))
p1
ggsave(p1, file="03_simulation_study/04_graphical_representation_prior/output_images/01_locations_neurons_nocl.pdf", width = 3, height=3)
ggsave(p1, file="03_simulation_study/04_graphical_representation_prior/output_images/01_locations_neurons_nocl.png", width = 3, height=3)



#---------# #-----------# #---------# #---------# #-----------# #---------#
#---------# #-----------#           PSBM          #-----------# #---------#
#---------# #-----------# #---------# #---------# #-----------# #---------#

couples = expand.grid(1:22, 1:22)

out = apply(couples, 1, function(x) sum((locations[x[1],] - locations[x[2],])^2)  )
df = data.frame(couples, out)
str(df)

library(paletteer)
p11 = ggplot(df, aes(factor(Var1), factor(Var2))) +
  geom_tile(aes(fill = out)) + 
  theme_minimal() +
  scale_fill_paletteer_c("grDevices::YlGnBu", direction=1)+
  theme(
    panel.background = element_rect(fill='transparent', color=NA), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    # panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    panel.border = element_rect(color = "darkgray ", fill=NA),
    # axis.line.y.left = element_line(color="gray"),
    axis.line.x.bottom = element_line(color="grey80"),
    #
    legend.position = "none",
    legend.background = element_rect(fill='transparent', color = 'transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent', color = 'transparent'), #transparent legend panel
    # legend.title=element_blank(),
    legend.text = element_text(size=12),
    strip.text = element_text(size=12),
    text = element_text(size = 12),
    axis.title.x=element_text(size=12),
    axis.title.y=element_text(size=12),
    strip.background = element_rect( fill=NA, color="grey" )
  ) +
  xlab("Neuron") + ylab("Neuron") +
  scale_y_discrete(limits=rev) +
  scale_x_discrete(breaks = seq(0,22,by=4)) +
  scale_y_discrete(breaks = seq(0,22,by=4)) 
p11
ggsave(p11, file="03_simulation_study/04_graphical_representation_prior/output_images/01_omega.pdf", width = 3, height=3)
ggsave(p11, file="03_simulation_study/04_graphical_representation_prior/output_images/01_omega.png", width = 3, height=3)


#---------# #-----------# #---------# #---------# #-----------# #---------#
#---------# #-----------#     latent processes    #-----------# #---------#
#---------# #-----------# #---------# #---------# #-----------# #---------#

step_fun = matrix(-2, 3, TT)
step_fun[1, c(1:4)] = -5
step_fun[1, c(20:23)] = 2
step_fun[1, c(24:25)] = 5
step_fun[1, c(49:53)] = 4

step_fun[2, 35] = 3
step_fun[2, 36:37] = 9
step_fun[2, 38] = 1

step_fun[3, 62:63] = 6
step_fun[3, 65:66] = 5

meanGP = t(apply(step_fun, 1, function(x) lm(x~poly(1:TT,7))$fitted.values ))*4

df = data.frame("time" = rep(1:TT, 3), 
                "p" = c(meanGP[1,], meanGP[2,], meanGP[3,]) ,
                "cl" = factor(rep(c("Atom 1","Atom 2","Atom 3"), each=TT), 
                              levels = c("Atom 3","Atom 2","Atom 1"))
)

p2 = ggplot(df, aes(x=time, y=p, color=cl)) +
  geom_line(lwd=1.2) + 
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
    legend.position = "none",
    legend.background = element_rect(fill='transparent', color = 'transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent', color = 'transparent'), #transparent legend panel
    # legend.title=element_blank(),
    legend.text = element_text(size=12),
    strip.text = element_text(size=12),
    text = element_text(size = 12),
    axis.title.x=element_text(size=12),
    axis.title.y=element_text(size=12),
    strip.background = element_rect( fill=NA, color="gray" )
  ) + 
  xlab("Time") + 
  ylab("Value") + 
  scale_color_manual(values =  c("deeppink3", "#00BFC4", "blue"))+
  facet_wrap(~cl, nrow=3) 
p2
ggsave(p2, file="03_simulation_study/04_graphical_representation_prior/output_images/02_GP.pdf", width = 3.5, height=4)
ggsave(p2, file="03_simulation_study/04_graphical_representation_prior/output_images/02_GP.png", width = 3.5, height=4)


df = data.frame("time" = rep(1:TT, 3), 
                "p" = pnorm(c(meanGP[1,], meanGP[2,], meanGP[3,])) ,
                "cl" = factor(rep(c("Atom 1","Atom 2","Atom 3"), each=TT), 
                              levels = c("Atom 3","Atom 2","Atom 1"))
)

p3 = ggplot(df, aes(x=time, y=p, color=cl)) +
  geom_line(lwd=1.2) + 
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
    legend.position = "none",
    legend.background = element_rect(fill='transparent', color = 'transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent', color = 'transparent'), #transparent legend panel
    # legend.title=element_blank(),
    legend.text = element_text(size=12),
    strip.text = element_text(size=12),
    text = element_text(size = 12),
    axis.title.x=element_text(size=12),
    axis.title.y=element_text(size=12),
    strip.background = element_rect( fill=NA, color="gray" )
  ) + 
  xlab("Time") + 
  ylab("Probability") + 
  scale_color_manual(values =  c("deeppink3", "#00BFC4", "blue"))+
  facet_wrap(~cl, nrow=3) 
ggsave(p3, file="03_simulation_study/04_graphical_representation_prior/output_images/03_probit_GP.pdf", width = 3.5, height=4)
ggsave(p3, file="03_simulation_study/04_graphical_representation_prior/output_images/03_probit_GP.png", width = 3.5, height=4)






#---------# #-----------# #---------# #---------# #-----------# #---------#
#---------# #-----------#          signal         #-----------# #---------#
#---------# #-----------# #---------# #---------# #-----------# #---------#

par_tau2 = 5
par_sigma2 = 1

# generate signal
set.seed(143)
sig1 = gen_signal(n, TT, cluster_neurons, meanGP = meanGP, par_tau2=par_tau2, par_sigma2 = par_sigma2)
signal = sig1$signal

library(tidyr)
long <- data.frame("neu" = 1:n, "t" = signal) %>% 
  pivot_longer(
    cols = !neu,
    names_to = "time",
    values_to = "value"
  )
long$t = rep(1:TT,n)
long$cl = rep(cluster_neurons, each = TT)
long$value = long$value*long$cl
str(long)

p4 = ggplot(long, aes(factor(t), factor(neu))) +
  geom_tile(aes(fill = factor(value)), colour = "grey80") + 
  scale_fill_manual(values = c("grey95", "blue", "#00BFC4", "deeppink3")) +
  scale_x_discrete(breaks = seq(0,80,by=20)) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill='transparent', color=NA), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    # panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    panel.border = element_rect(color = "darkgray ", fill=NA),
    # axis.line.y.left = element_line(color="gray"),
    axis.line.x.bottom = element_line(color="grey80"),
    #
    legend.position = "none",
    legend.background = element_rect(fill='transparent', color = 'transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent', color = 'transparent'), #transparent legend panel
    # legend.title=element_blank(),
    legend.text = element_text(size=12),
    strip.text = element_text(size=12),
    text = element_text(size = 12),
    axis.title.x=element_text(size=12),
    axis.title.y=element_text(size=12),
    strip.background = element_rect( fill=NA, color="grey" )
  ) +
  xlab("Time") + ylab("Neuron")
p4
ggsave(p4, file="03_simulation_study/04_graphical_representation_prior/output_images/04_activations.pdf", width = 4, height=5)
ggsave(p4, file="03_simulation_study/04_graphical_representation_prior/output_images/04_activations.png", width = 4, height=5)






#---------# #-----------# #---------# #---------# #-----------# #---------#
#---------# #-----------#       amplitudes        #-----------# #---------#
#---------# #-----------# #---------# #---------# #-----------# #---------#

w0 = 0.5

set.seed(123)
atoms = rgamma(50, 5, 0.5)
atoms

idx= sort(atoms[1:7], index.return=T, decreasing=T)$ix

sticks = rbeta(50, 1, 1)
w = numeric(50)
w[1] = sticks[1]
w[2:50] = sapply(2:50, function(x) sticks[x] * prod(1-sticks[1:(x-1)]) )


df <- data.frame(x = c(0,atoms[idx]), 
                 y = c(w0,((1-w0)*w[1:7])[idx]),
                 col= factor(1:8))

# plot
p5 = ggplot(df, aes(x=x, y=y)) +
  geom_segment( aes(x=x, xend=x, y=0, yend=y, color=col),lwd=0.7) +
  geom_point( size=2, aes(fill=col, color=col), alpha=0.7, shape=20, stroke=1) +
  scale_color_manual(values = c("grey55", inferno(11)[9:3])) +
  scale_fill_manual(values = c("grey55", inferno(11)[9:3])) +
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
    legend.position = "none",
    legend.background = element_rect(fill='transparent', color = 'transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent', color = 'transparent'), #transparent legend panel
    # legend.title=element_blank(),
    legend.text = element_text(size=12),
    strip.text = element_text(size=12),
    text = element_text(size = 12),
    axis.title.x=element_text(size=12),
    axis.title.y=element_text(size=12),
    strip.background = element_rect( fill=NA, color="gray" )
  ) +
  xlab("Amplitude") + ylab("Probability") +
  xlim(0,20) + 
  ylim(0,0.55) 

p5
ggsave(p5, file="03_simulation_study/04_graphical_representation_prior/output_images/05_distribution_ampl.pdf", width = 3, height=3)
ggsave(p5, file="03_simulation_study/04_graphical_representation_prior/output_images/05_distribution_ampl.png", width = 3, height=3)



set.seed(123)
unique_amplitudes = atoms[sample(1:50, prob=w, replace=T)]



# # spike amplitudes
# unique_amplitudes = c(5,
#                       6,6,6,
#                       10,10,10,10,
#                       14,
#                       25) 

# given the signal, generate the calcium traces 
cal = gen_amplitudes_trace(signal, ampl = unique_amplitudes, gamma = 0.9, tau2 = 1, sigma2 = 3) 

sim = list(n = n, TT = TT,
           calcium = cal$trace,
           loc = locations,
           cluster_neurons = cluster_neurons,
           signal = signal,
           amplitudes = cal$amplitudes,
           cluster_a = cal$clus_ampl,
           unique_amplitudes = unique_amplitudes,
           latent_prob = atoms,
           atoms_GP = sig1$atoms_GP
)






#---------# #-----------# #---------# #---------# #-----------# #---------#
#---------# #-----------#   signal + amplitudes   #-----------# #---------#
#---------# #-----------# #---------# #---------# #-----------# #---------#



long <- data.frame("neu" = 1:n, "t" = sim$amplitudes) %>% 
  pivot_longer(
    cols = !neu,
    names_to = "time",
    values_to = "value"
  )
long$t = rep(1:TT,n)
table(long$value)

library(viridis)
p6 = ggplot(long, aes(factor(t), factor(neu))) +
  geom_tile(aes(fill = factor(value)), colour = "grey80") + 
  scale_x_discrete(breaks = seq(0,80,by=20)) +
  scale_fill_manual(values = c("grey95", inferno(11)[9:3])) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill='transparent', color=NA), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    # panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    panel.border = element_rect(color = "darkgray ", fill=NA),
    # axis.line.y.left = element_line(color="gray"),
    axis.line.x.bottom = element_line(color="grey80"),
    #
    legend.position = "none",
    legend.background = element_rect(fill='transparent', color = 'transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent', color = 'transparent'), #transparent legend panel
    # legend.title=element_blank(),
    legend.text = element_text(size=12),
    strip.text = element_text(size=12),
    text = element_text(size = 12),
    axis.title.x=element_text(size=12),
    axis.title.y=element_text(size=12),
    strip.background = element_rect( fill=NA, color="gray" )
  ) +
  xlab("Time") + ylab("Neuron")

p6
ggsave(p6, file="03_simulation_study/04_graphical_representation_prior/output_images/06_spike_ampl.pdf", width = 4, height=5)
ggsave(p6, file="03_simulation_study/04_graphical_representation_prior/output_images/06_spike_ampl.png", width = 4, height=5)





#---------# #-----------# #---------# #---------# #-----------# #---------#
#---------# #-----------#          calcium        #-----------# #---------#
#---------# #-----------# #---------# #---------# #-----------# #---------#

cal = sim$calcium
cal = cal/( max(cal)-30)

calcium_active = t(cal) + matrix(rep(1:n, TT), TT, n, byrow=T)
calcium_active = reshape2::melt((calcium_active))
calcium_active$time = rep(1:TT, n)
calcium_active$neuron = as.factor(as.numeric(calcium_active$Var2))
calcium_active$spike = as.numeric(reshape2::melt(t(sim$signal))$value)
calcium_active$Cluster = as.factor(rep(cluster_neurons, each = TT))
calcium_active = rev(calcium_active)
calcium_active$rev_neuron = rev(calcium_active$neuron)

p7 = ggplot(data = calcium_active, aes(x = time, y = value, color = Cluster)) +
  geom_ribbon(aes(x = time, ymax = value, ymin =as.numeric(neuron), fill=Cluster, color=Cluster, group = rev_neuron), alpha=0.08, lwd=0.1,
              show.legend=TRUE)+
  geom_line(aes(x = time, y = value, color=Cluster, group = rev_neuron), alpha=1, lwd=0.3,
            show.legend=TRUE)  +
  geom_point(data = calcium_active[calcium_active$spike>0,], aes(x = time, y = as.numeric(neuron)), alpha = 0.35, size = 0.9,
             show.legend=TRUE) +
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
    legend.text = element_text(size=10),
    strip.text = element_text(size=10),
    strip.background = element_rect( fill=NA, color="gray" )
  )+
  guides(colour = guide_legend(nrow = 1)) +
  xlab("Time")  +
  ylab("Neuron") +
  scale_fill_manual(values = c("blue", "#00BFC4", "deeppink3")) +
  scale_color_manual(values = c("blue", "#00BFC4", "deeppink3"))

p7
ggsave(p7, file="03_simulation_study/04_graphical_representation_prior/output_images/07_calcium.pdf", width = 4.5, height=6)
ggsave(p7, file="03_simulation_study/04_graphical_representation_prior/output_images/07_calcium.png", width = 4.5, height=6)

