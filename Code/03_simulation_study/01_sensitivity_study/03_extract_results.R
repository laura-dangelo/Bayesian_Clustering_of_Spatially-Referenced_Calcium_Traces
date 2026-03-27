
#-------------------# #-------------------# #-------------------#
##                EXTRACT QUANTITIES OF INTEREST               ##
#-------------------# #-------------------# #-------------------#

# These scripts replicate the sensitivity study reported in Section S4.2 of the Supplementary Material.
# Specifically, this script takes as input the runs of the Gibbs sampler and compute the classification error rates and ARI.
# The outputs are saved in the output_RDS folder.

#  _________________________________
#  YOU CAN AVOID RUNNING THIS SCRIPT 
#  ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
# the RDS files produced by this script are available in the "Code/03_simulation_study/01_sensitivity_study/output_RDS" folder.


library(salso)
library(mclust)

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
    nameopen = paste0("03_simulation_study/01_sensitivity_study/simulated_data/simulated_data", n_sim, ".RDS")
    sim = readRDS(nameopen)

    name = paste0("03_simulation_study/01_sensitivity_study/results/run_gibbs_alow_", 
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


saveRDS(est_sigma2, file = "03_simulation_study/01_sensitivity_study/output_RDS/est_sigma2.RDS")
saveRDS(est_tau2, file = "03_simulation_study/01_sensitivity_study/output_RDS/est_tau2.RDS")
saveRDS(errors_spike_detection, file = "03_simulation_study/01_sensitivity_study/output_RDS/errors_spike_detection.RDS")
saveRDS(false_positive_rates, file = "03_simulation_study/01_sensitivity_study/output_RDS/false_positive_rates.RDS")
saveRDS(false_negative_rates, file = "03_simulation_study/01_sensitivity_study/output_RDS/false_negative_rates.RDS")
saveRDS(aris_neuron, file = "03_simulation_study/01_sensitivity_study/output_RDS/ARIs.RDS")


est_sigma2 = readRDS( file = "03_simulation_study/01_sensitivity_study/output_RDS/est_sigma2.RDS")
est_tau2= readRDS( file = "03_simulation_study/01_sensitivity_study/output_RDS/est_tau2.RDS")
errors_spike_detection= readRDS( file = "03_simulation_study/01_sensitivity_study/output_RDS/errors_spike_detection.RDS")
false_positive_rates= readRDS( file = "03_simulation_study/01_sensitivity_study/output_RDS/false_positive_rates.RDS")
false_negative_rates= readRDS( file = "03_simulation_study/01_sensitivity_study/output_RDS/false_negative_rates.RDS")
aris_neuron= readRDS( file = "03_simulation_study/01_sensitivity_study/output_RDS/ARIs.RDS")
