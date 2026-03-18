
library(salso)
library(ggplot2)
library(viridisLite)
library(mclust)

rm(list=ls())
replications_sim = 50

seq_alow = c(0, 0.5, 1)
id_pars_gamma = c(1,2,3)
pars_gamma = matrix(c(3, 0.1, 4, 1, 10, 1), 3, 2, byrow = T)

parameter_comb = expand.grid(x = seq_alow, y = id_pars_gamma)
parameter_comb[,3] = sapply(parameter_comb[,2], function(x) pars_gamma[x,2])
parameter_comb[,2] = sapply(parameter_comb[,2], function(x) pars_gamma[x,1])
parameter_comb


compute_class_error = function(true_sig, estimate_sig)
{
  diff_error = true_sig - estimate_sig
  # 0 : correct
  # 1 : false negative
  # -1 : false positive
  
  tot = nrow(true_sig) * ncol(true_sig)
  false_positives = sum(diff_error == -1)
  false_negatives = sum(diff_error == 1)
  tot_positives = sum(true_sig == 1)
  tot_negatives = sum(true_sig == 0)
  
  FPR = false_positives / tot_negatives
  FNR = false_negatives / tot_positives
  
  class_error = (false_positives + false_negatives) / tot
  
  list("FPR" = FPR, "FNR" = FNR, "class_error" = class_error, "tot_FP" = false_positives, "tot_FN" = false_negatives)
}



est_sigma2 = matrix(0, replications_sim, nrow(parameter_comb))
est_tau2 =  matrix(0, replications_sim, nrow(parameter_comb))
aris_neuron =  matrix(0, replications_sim, nrow(parameter_comb))

errors_spike_detection = matrix(0, replications_sim, nrow(parameter_comb))
false_positive_rates = matrix(0, replications_sim, nrow(parameter_comb))
false_negative_rates = matrix(0, replications_sim, nrow(parameter_comb))

colnames(est_sigma2) = colnames(est_tau2) = 
colnames(aris_neuron) = 
colnames(errors_spike_detection) = 
colnames(false_positive_rates) = 
colnames(false_negative_rates) = 
c("a_low = 0, Gamma(3, 0.1)", "a_low = 0.5, Gamma(3, 0.1)", "a_low = 1, Gamma(3, 0.1)", 
  "a_low = 0, Gamma(4, 1)", "a_low = 0.5, Gamma(4, 1)", "a_low = 1, Gamma(4, 1)",  
  "a_low = 0, Gamma(10, 1)", "a_low = 0.5, Gamma(10, 1)", "a_low = 1, Gamma(10, 1)" )

for(p in 1:nrow(parameter_comb)){
  for(n_sim in 1:replications_sim){
    nameopen = paste0("10_simulations_sensitivity/simulated_data/simulated_data", n_sim, ".RDS")
    sim = readRDS(nameopen)

    name = paste0("10_simulations_sensitivity/results/run_gibbs_alow_", 
                  parameter_comb[p,1], "_", parameter_comb[p,2], "_", parameter_comb[p,3], 
                  "_simu", n_sim, ".RDS")
    out = readRDS(name)
    
    TT = dim(out$latent_signal)[1]
    N = dim(out$AA)[2]
    
    est_sigma2[n_sim, p] = mean(out$sigma2)
    est_tau2[n_sim, p] = mean(out$tau2)
    
    est_cluster_neurons = salso(t(out$cluster_signal+1) , maxZealousAttempts = 100,
                                maxNClusters = 25,
                                nCores = 20,
                                nRuns = 500)
    aris_neuron[n_sim, p] =  adjustedRandIndex(est_cluster_neurons, sim$cluster_neurons)

    
    PPS = matrix(0, TT, N)
    for(i in 1:N) {
      PPS[,i] = apply(out$AA[,i,1:length(out$gamma)], 1, function(x) mean(x>0))
    }
    threshold = 0.8
    
    estimated_spikes = (PPS>threshold)
    err = compute_class_error(t(sim$amplitudes>0), estimated_spikes)
    
    errors_spike_detection[n_sim, p] = err$class_error
    false_positive_rates[n_sim, p] = err$FPR
    false_negative_rates[n_sim, p] = err$FNR
    cat(paste("nsim:",n_sim,"par_comb:", p,"\n"))
    
  }
}

est_sigma2
est_tau2
errors_spike_detection
false_positive_rates
false_negative_rates
aris_neuron


saveRDS(est_sigma2, file = "10_simulations_sensitivity/output_RDS/est_sigma2.RDS")
saveRDS(est_tau2, file = "10_simulations_sensitivity/output_RDS/est_tau2.RDS")
saveRDS(errors_spike_detection, file = "10_simulations_sensitivity/output_RDS/errors_spike_detection.RDS")
saveRDS(false_positive_rates, file = "10_simulations_sensitivity/output_RDS/false_positive_rates.RDS")
saveRDS(false_negative_rates, file = "10_simulations_sensitivity/output_RDS/false_negative_rates.RDS")
saveRDS(aris_neuron, file = "10_simulations_sensitivity/output_RDS/ARIs.RDS")


est_sigma2 = readRDS( file = "10_simulations_sensitivity/output_RDS/est_sigma2.RDS")
est_tau2= readRDS( file = "10_simulations_sensitivity/output_RDS/est_tau2.RDS")
errors_spike_detection= readRDS( file = "10_simulations_sensitivity/output_RDS/errors_spike_detection.RDS")
false_positive_rates= readRDS( file = "10_simulations_sensitivity/output_RDS/false_positive_rates.RDS")
false_negative_rates= readRDS( file = "10_simulations_sensitivity/output_RDS/false_negative_rates.RDS")
aris_neuron= readRDS( file = "10_simulations_sensitivity/output_RDS/ARIs.RDS")




#----------# #----------# #----------# #----------#
#----------# PLOT ARI

alow_seq = factor(rep(c(0,0.5,1), replications_sim))
params_seq = factor(rep(c("a301", "b41", "c101"), each = replications_sim))

dfari = data.frame("ARI" = c(aris_neuron), "a_low" = alow_seq, "gamma_params" = params_seq)


ggplot(dfari, aes(y = ARI, x = a_low)) +
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
  facet_wrap(gamma_params ~ ., 
             labeller = labeller(gamma_params = c("a301" = "Gamma(3, 0.1)",
                                                   "b41" = "Gamma(4, 1)",
                                                   "c101" = "Gamma(10, 1)") ) )+
  xlab(expression(bar(a)))  +
  ylab("Adjusted Rand index")


ggsave("10_simulations_sensitivity/output_images/adjustedRandIndex.pdf", width = 7, height = 4)
ggsave("10_simulations_sensitivity/output_images/adjustedRandIndex.eps", width = 7, height = 4, device = cairo_ps)




#----------# #----------# #----------# #----------#
#----------# PLOT ERROR RATE OF SPIKE DETECTION

#----------#  Misclassification error

alow_seq = factor(rep(c(0,0.5,1), replications_sim))
params_seq = factor(rep(c("a301", "b41", "c101"), each = replications_sim))

dferror = data.frame("error" = c(errors_spike_detection), "a_low" = alow_seq, "gamma_params" = params_seq)


ggplot(dferror, aes(y = error, x = a_low, fill=a_low)) +
  geom_boxplot(alpha = 0.3) +
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
  facet_wrap(gamma_params ~ ., 
             labeller = labeller(gamma_params = c("a301" = "Gamma(3, 0.1)",
                                                  "b41" = "Gamma(4, 1)",
                                                  "c101" = "Gamma(10, 1)") ) )+
  xlab(expression(paste(bar(a))))  +
  ylab("Misclassification error rate")

ggsave("10_simulations_sensitivity/output_images/error_rates.pdf", width = 7, height = 4)
ggsave("10_simulations_sensitivity/output_images/error_rates.eps", width = 7, height = 4, device = cairo_ps)





#----------#  False negative rate
alow_seq = factor(rep(c(0,0.5,1), replications_sim))
params_seq = factor(rep(c("a301", "b41", "c101"), each = replications_sim))

dffn = data.frame("error" = c(false_negative_rates), "a_low" = alow_seq, "gamma_params" = params_seq)


ggplot(dffn, aes(y = error, x = a_low, fill=a_low)) +
  geom_boxplot(alpha = 0.3) +
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
  facet_wrap(gamma_params ~ ., 
             labeller = labeller(gamma_params = c("a301" = "Gamma(3, 0.1)",
                                                  "b41" = "Gamma(4, 1)",
                                                  "c101" = "Gamma(10, 1)") ) )+
  xlab(expression(paste(bar(a))))  +
  ylab("False negative rate")

ggsave("10_simulations_sensitivity/output_images/fn_rates.pdf", width = 5, height = 4)
ggsave("10_simulations_sensitivity/output_images/fn_rates.eps", width = 5, height = 4, device = cairo_ps)



#----------#  False positive rate

alow_seq = factor(rep(c(0,0.5,1), replications_sim))
params_seq = factor(rep(c("a301", "b41", "c101"), each = replications_sim))

dffp = data.frame("error" = c(false_positive_rates), "a_low" = alow_seq, "gamma_params" = params_seq)


ggplot(dffn, aes(y = error, x = a_low, fill=a_low)) +
  geom_boxplot(alpha = 0.3) +
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
  facet_wrap(gamma_params ~ ., 
             labeller = labeller(gamma_params = c("a301" = "Gamma(3, 0.1)",
                                                  "b41" = "Gamma(4, 1)",
                                                  "c101" = "Gamma(10, 1)") ) )+
  xlab(expression(paste(bar(a))))  +
  ylab("False positive rate")

ggsave("10_simulations_sensitivity/output_images/fp_rates.pdf", width = 5, height = 4)
ggsave("10_simulations_sensitivity/output_images/fp_rates.eps", width = 5, height = 4, device = cairo_ps)



#----------#  together
alow_seq = factor(rep(c(0,0.5,1), replications_sim))
params_seq = factor(rep(c("a301", "b41", "c101"), each = replications_sim))

errors_full = data.frame(rbind(dffp, dffn, dferror))
colnames(errors_full)[1] = "value"
errors_full$Error = c(cbind(rep("False positive", nrow(dffp)),rep("False negative", nrow(dffp)),rep("Misclassification", nrow(dffp)) ))


ggplot(errors_full, aes(y = value, x = a_low)) +
  geom_boxplot(aes(fill=Error),  alpha = 0.3) +
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
  facet_wrap(gamma_params ~ ., 
             labeller = labeller(gamma_params = c("a301" = "Gamma(3, 0.1)",
                                                  "b41" = "Gamma(4, 1)",
                                                  "c101" = "Gamma(10, 1)") ) )+
  xlab(expression(paste(bar(a))))  +
  ylab("Error rate")

ggsave("10_simulations_sensitivity/output_images/rates_complete.pdf", width = 7, height = 4)
ggsave("10_simulations_sensitivity/output_images/rates_complete.eps", width = 7, height = 4, device = cairo_ps)

ggsave("10_simulations_sensitivity/output_images/rates_complete_wide.pdf", width = 10, height = 4)
ggsave("10_simulations_sensitivity/output_images/rates_complete_wide.eps", width = 10, height = 4, device = cairo_ps)

