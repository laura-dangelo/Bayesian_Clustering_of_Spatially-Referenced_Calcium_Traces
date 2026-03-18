
#---------# #---------# #---------# #---------# #---------# 
#---------#     Toy example from probit-SB      #---------# 
#---------# #---------# #---------# #---------# #---------# 

library(tidyverse)
library(mvtnorm)
library(viridisLite)

trunc <- 50


#---------# #---------# #---------# #---------# 
#---------#     Generate data       #---------# 
set.seed(123)
L <- replicate(2,rnorm(100,1,0.5))
L[1:30,] <- L[1:30,]+3
L[1:30+30,] <- L[1:30+30,]-3
truecl <- rep(1:3,c(30,30,40))

colneu = rep(1, nrow(L))
colneu[c(1,4,31,36,61,67)] = 2:7

sizeneu = rep(2, nrow(L))
sizeneu[c(1,4,31,36,61,67)] = 2.7

df = data.frame(cbind(L,truecl,factor(colneu),sizeneu))

colors = c("black", "firebrick", "red", "forestgreen", "darkseagreen", "turquoise4", "blue")
# locations
ggplot(data=df, aes(x=V1, y=V2, col=as.factor(V4))) + 
  geom_point(size = sizeneu, alpha = 0.7) +
  scale_color_manual(name="cluster", values = colors)+
  theme_minimal() +
  theme(
    panel.background = element_rect(fill='transparent', color=NA), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    # panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    panel.border = element_rect(color = "darkgray ", fill=NA),
    # axis.line.y.left = element_line(color="gray"),
    axis.line.x.bottom = element_line(color="gray"),
    axis.text = element_text(size=11),
    axis.title = element_text(size=13),
    #
    legend.position = "none",
    legend.background = element_rect(fill='transparent', color = 'transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent', color="transparent"), #transparent legend panel
    legend.text = element_text(size=10),
    strip.text = element_text(size=13),
    strip.background = element_rect( fill=NA, color="gray" )
  )+
  xlab("x coordinate")  +
  ylab("y coordinate") 


ggsave("02_data_analysis/03_bSCDC_sensitivity_study/output_images/01_points_space.pdf", width = 6, height = 4.5)


#---------# #---------# #---------# #---------# #---------# #---------# #---------# #---------# 

#---------# #---------# #---------# #---------# #---------# 
#---------#     Probit-SB assuming theta = 1    #---------# 

theta <- 1
tau <- 1
cov <- tau * exp(-.5 * as.matrix(dist(L)^2) * theta) 
# covariance matrix
pheatmap::pheatmap(cov, cluster_rows = F,cluster_cols = F, #fontsize_row=0.2,fontsize_col=0.2,
                   labels_row= "j", labels_col= "i", angle_col = "0")


# generate 30 replicates of n-variates normal with this covariance matrix
set.seed(1234)
alpha <- rmvnorm(trunc,
                 rep(0,100),
                 sigma = cov)
matplot(t(alpha)[,1:4], type="l", xlab = "x", ylab = "y", lty=1)

pheatmap::pheatmap((alpha),cluster_rows = F,cluster_cols = F)


psbm <- function(a){
  n <- length(a)
  v <- pnorm(a)
  one_v <- 1-v
  W <- v
  W[2:n] <- v[2:n]*cumprod(one_v[1:(n-1)])
  return(W)
}


# 1,4,31,36,61,67
orderneu = paste("Neuron", sort(rep(c(1,4,31,36,61,67), 15)))
orderneu = factor(orderneu, levels = c("Neuron 1", "Neuron 31", "Neuron 61",
                                    "Neuron 4", "Neuron 36", "Neuron 67"))               

colorneu2 = sort(rep(1:6, 15)  )               

colors = c("firebrick", "red", "springgreen4", "chartreuse3", "turquoise3", "blue")

ggplot(data.frame("x"=rep(1:15,3), 
                  "prob"= c( psbm(alpha[,1])[1:15],
                             psbm(alpha[,4])[1:15],
                             psbm(alpha[,31])[1:15],
                             psbm(alpha[,36])[1:15],
                             psbm(alpha[,61])[1:15],
                             psbm(alpha[,67])[1:15]
                             ), 
                  "id" =  orderneu ,
                  "col" = factor(colorneu2)
                  ), aes(x=x, y=prob)) + 
  geom_point(aes(col=col), alpha=.5) + 
  geom_segment(aes(x=x, xend=x, y=0, yend=prob, col=col)) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill='transparent', color=NA), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    # panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    panel.border = element_rect(color = "darkgray ", fill=NA),
    # axis.line.y.left = element_line(color="gray"),
    axis.line.x.bottom = element_line(color="gray"),
    axis.text = element_text(size=8),
    axis.title = element_text(size=11),
    #
    legend.position = "none",
    legend.background = element_rect(fill='transparent', color = 'transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent', color="transparent"), #transparent legend panel
    legend.text = element_text(size=10),
    strip.text = element_text(size=10),
    strip.background = element_rect( fill=NA, color="gray" )
  )+
  xlab("i")  +
  ylab("Probability") +
  scale_colour_manual(values=colors) + 
  scale_x_continuous(breaks=1:15, labels = ) +
  facet_wrap(~id, ncol=3)

ggsave("02_data_analysis/03_bSCDC_sensitivity_study/output_images/04_prob_theta1.pdf", width = 9, height = 5)





#---------# #---------# #---------# #---------# #---------# 
#---------#   Probit-SB assuming theta = 1000   #---------# 

theta <- 1000
tau <- 1
cov <- tau * exp(-.5 * as.matrix(dist(L)^2)*theta) 
# covariance matrix
pheatmap::pheatmap(cov, cluster_rows = F,cluster_cols = F, #fontsize_row=0.2,fontsize_col=0.2,
                   labels_row= "j", labels_col= "i", angle_col = "0")


# generate 30 replicates of n-variates normal with this covariance matrix
set.seed(12)
alpha <- rmvnorm(trunc,
                 rep(0,100),
                 sigma = cov)
matplot(t(alpha)[,1:3], type="l", xlab = "x", ylab = "y", lty=1)

pheatmap::pheatmap((alpha),cluster_rows = F,cluster_cols = F)


psbm <- function(a){
  n <- length(a)
  v <- pnorm(a)
  one_v <- 1-v
  W <- v
  W[2:n] <- v[2:n]*cumprod(one_v[1:(n-1)])
  return(W)
}

# 1,4,31,36,61,67
orderneu = paste("Neuron", sort(rep(c(1,4,31,36,61,67), 15)))
orderneu = factor(orderneu, levels = c("Neuron 1", "Neuron 31", "Neuron 61",
                                       "Neuron 4", "Neuron 36", "Neuron 67"))               

colorneu2 = sort(rep(1:6, 15)  )               

colors = c("firebrick", "red", "springgreen4", "chartreuse3", "turquoise3", "blue")

ggplot(data.frame("x"=rep(1:15,3), 
                  "prob"= c( psbm(alpha[,1])[1:15],
                             psbm(alpha[,4])[1:15],
                             psbm(alpha[,31])[1:15],
                             psbm(alpha[,36])[1:15],
                             psbm(alpha[,61])[1:15],
                             psbm(alpha[,67])[1:15]
                  ), 
                  "id" =  orderneu ,
                  "col" = factor(colorneu2)
), aes(x=x, y=prob)) + 
  geom_point(aes(col=col), alpha=.5) + 
  geom_segment(aes(x=x, xend=x, y=0, yend=prob, col=col)) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill='transparent', color=NA), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    # panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    panel.border = element_rect(color = "darkgray ", fill=NA),
    # axis.line.y.left = element_line(color="gray"),
    axis.line.x.bottom = element_line(color="gray"),
    axis.text = element_text(size=8),
    axis.title = element_text(size=11),
    #
    legend.position = "none",
    legend.background = element_rect(fill='transparent', color = 'transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent', color="transparent"), #transparent legend panel
    legend.text = element_text(size=10),
    strip.text = element_text(size=10),
    strip.background = element_rect( fill=NA, color="gray" )
  )+
  xlab("i")  +
  ylab("Probability") +
  scale_colour_manual(values=colors) + 
  scale_x_continuous(breaks=1:15, labels = ) +
  facet_wrap(~id, ncol=3)

ggsave("02_data_analysis/03_bSCDC_sensitivity_study/output_images/05_prob_theta1000.pdf", width = 9, height = 5)


