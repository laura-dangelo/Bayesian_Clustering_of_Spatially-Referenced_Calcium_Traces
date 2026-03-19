
#-------------------# #-------------------# #-------------------#
##                EXTRACT QUANTITIES OF INTEREST               ##
#-------------------# #-------------------# #-------------------#

# These scripts replicate the comparison with the two-step approach reported in Section S4.3 of the Supplementary Material.
# Specifically, this script computes the classification error rates.

library(salso)
library(mclust)
replications_sim = 50


compute_class_error_JW = function(true_sig, estimated_signal)
{
  diff_error = true_sig - estimated_signal
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


errors_spike_detection_JW = rep(0, replications_sim)
false_positive_rates_JW = rep(0, replications_sim)
false_negative_rates_JW = rep(0, replications_sim)


for(n_sim in 1:replications_sim){
    nameopen = paste0("03_simulation_study/01_sensitivity_study/simulated_data/simulated_data", n_sim, ".RDS")
    sim = readRDS(nameopen)

    name = paste0("03_simulation_study/02_comparison_two-step_synthetic_data/results/estimates_JW_simu", n_sim, ".RDS")
    out = readRDS(file = name)
    
    TT = ncol(out)
    N = nrow(out)
    
    err = compute_class_error_JW((sim$amplitudes>0), out)
    
    errors_spike_detection_JW[n_sim] = err$class_error
    false_positive_rates_JW[n_sim] = err$FPR
    false_negative_rates_JW[n_sim] = err$FNR
    
}

errors_spike_detection_JW
false_positive_rates_JW
false_negative_rates_JW


saveRDS(errors_spike_detection_JW, file = "03_simulation_study/02_comparison_two-step_synthetic_data/output_RDS/errors_spike_detection_JW.RDS")
saveRDS(false_positive_rates_JW, file = "03_simulation_study/02_comparison_two-step_synthetic_data/output_RDS/false_positive_rates_JW.RDS")
saveRDS(false_negative_rates_JW, file = "03_simulation_study/02_comparison_two-step_synthetic_data/output_RDS/false_negative_rates_JW.RDS")

