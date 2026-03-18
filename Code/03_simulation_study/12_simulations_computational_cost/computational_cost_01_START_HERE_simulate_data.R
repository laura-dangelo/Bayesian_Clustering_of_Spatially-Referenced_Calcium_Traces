library(mvtnorm)
source("10_simulations_sensitivity/sensitivity00_functions_simulate_data.R")

#-------------------# #-------------------# #-------------------#
##                   GENERATE SIMULATED DATA
#-------------------# #-------------------# #-------------------#

# replications_sim: number of simulation replicates
replications_sim = 50

# n: number of neurons
# TT: length of the time series

seq_n = c(50,100,200,400) 
seq_T = c(100,200,400)

parameter_comb = expand.grid(x = seq_n, y = seq_T)
parameter_comb


for(i in 1:nrow(parameter_comb)){
  
  n = parameter_comb[i,1]
  TT = parameter_comb[i,2]
  
  for(n_sim in 1:replications_sim){
  
    set.seed(n_sim)
    
    locations = rbind(
      rmvnorm(round(n/4), c(0, 200), diag(2)*1000),
      rmvnorm(n - round(n/3)- round(n/4), c(110, 220), diag(2)*500),
      rmvnorm(round(n/3), c(80, 150), diag(2)*700)
    )
    cluster_neurons = c( rep(1, round(n/4)), rep(2, n - round(n/3)- round(n/4)), rep(3, round(n/3)) )
    index_mixed = sample(1:n, n*(25/100), replace=F)
    for(i in 1:length(index_mixed)){
      cluster_neurons[index_mixed[i]] = sample( setdiff(1:3, cluster_neurons[index_mixed[i]]), 1 )
    }
    # plot(locations, cex=2,  pch=21, bg=c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")[cluster_neurons])
    
    
    
    # parameters of the latent GP
    # length scale tau2: how wiggly are observations (amplitude w.r.t. t)
    # variance sigma2: amplitude of observations (amplitude w.r.t y)
    step_fun = matrix(-2, 3, TT)
    step_fun[1, c(1:4)] = -5
    int1 = round(TT/5):(round(TT/5)+5)
    step_fun[1, int1] = 3
    int2 = round(TT/2):(round(TT/2)+5)
    step_fun[1, int2] = 4
    
    int3 = round(TT/2):(round(TT/2)+1)
    step_fun[2, int3] = 3
    int4 = (round(TT/2)+2):(round(TT/2)+3)
    step_fun[2, int4] = 9
    int5 = (round(TT/2)+4):(round(TT/2)+6)
    step_fun[2, int5] = 1
    
    int6 = (round(3*TT/4)+1):(round(3*TT/4)+4)
    step_fun[3, int6] = 7
    int7 = (round(3*TT/4)+8):(round(3*TT/4)+9)
    step_fun[3, int5] = 2
    
    meanGP = t(apply(step_fun, 1, function(x) lm(x~poly(1:TT,7))$fitted.values ))*4
  
    
    # par(mfrow=c(3,1))
    # for(i in 1:3){
    #   plot(pnorm(meanGP[i,]),type="l")
    # }
    # par(mfrow=c(1,1))
    
    
    par_tau2 = 5
    par_sigma2 = 1
    
    # generate signal
    sig1 = gen_signal(n, TT, cluster_neurons, meanGP = meanGP, par_tau2=par_tau2, par_sigma2 = par_sigma2)
    signal = sig1$signal
    atoms = pnorm(sig1$atoms_GP)
    
    # spike smplitudes
    unique_amplitudes = c(5,
                          6,6,6,
                          10,10,10,10,
                          14,
                          25) 
    
    # given the signal, generate the calcium traces 
    cal = gen_amplitudes_trace(signal, ampl = unique_amplitudes, gamma = 0.9, tau2 = 1, sigma2 = 1.5) 
    
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
    name = paste0("12_simulations_computational_cost/simulated_data/simulated_data_n",n,"_TT", TT, "_sim", n_sim, ".RDS")
  
    saveRDS(sim, file = name)
  
  }
}
