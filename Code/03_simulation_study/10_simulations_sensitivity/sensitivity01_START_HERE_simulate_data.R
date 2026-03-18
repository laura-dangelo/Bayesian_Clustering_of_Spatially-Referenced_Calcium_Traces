library(mvtnorm)
source("10_simulations_sensitivity/sensitivity00_functions_simulate_data.R")

#-------------------# #-------------------# #-------------------#
##                   GENERATE SIMULATED DATA
#-------------------# #-------------------# #-------------------#

# replications_sim: number of simulation replicates
replications_sim = 50

# n: number of neurons
# TT: length of the time series
n = 100
TT = 80 # average length of the windows on real data is 78



for(n_sim in 1:replications_sim)
{

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
  step_fun[1, c(20:25)] = 3
  step_fun[1, c(49:53)] = 4
  
  step_fun[2, 45] = 3
  step_fun[2, 46:47] = 9
  step_fun[2, 48] = 1
  
  step_fun[3, 62:63] = 7
  step_fun[3, 65:66] = 2
  
  
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
  
  # par(mfrow=c(3,1))
  # for(i in 1:3){
  #   plot((atoms[i,]),type="l", ylim=c(0,1))  # evolution of the probability
  # }
  # par(mfrow=c(1,1))
  
  
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
  name = paste0("10_simulations_sensitivity/simulated_data/simulated_data", n_sim, ".RDS")

  saveRDS(sim, file = name)

}




#---------# #-----------# #---------# #---------# #-----------# #---------#
#---------# #-----------# #---------# #---------# #-----------# #---------#
##### example plot of simulated data 
# 
# n_sim=1
# nameopen = paste0("10_simulations_sensitivity/simulated_data/simulated_data", n_sim,".RDS")
# sim = readRDS(nameopen)
# 
# 
# 
# n = sim$n
# par(mfrow=c(1,2))
# const = 30
# nlines = n/2
# j=1
# plot(1:length(sim$calcium[1,]), sim$calcium[1,]/(const), type="n", ylim = c(1,nlines+1),
#      xlab = "Time", ylab = "sim$calcium level")
# idx = sort(sim$cluster_neurons, index.return = T)$ix
# j=1
# for(t in which(sim$amplitudes[idx[j],]>0)) { segments(t, j, t, j+0.2, col = sim$cluster_neurons[idx[j]] )}
# lines(sim$calcium[idx[j],]/(const)+j, col = sim$cluster_neurons[idx[j]])
# for(j in 2:(nlines)) {
#   for(t in which(sim$amplitudes[idx[j],]>0)) { segments(t, j, t, j+0.2, col = sim$cluster_neurons[idx[j]]) }
#   lines(sim$calcium[idx[j],]/(const)+j, col = sim$cluster_neurons[idx[j]])
# }
# plot(1:length(sim$calcium[j+1,]), sim$calcium[j+1,]/(const), type="n", ylim = c(1,nlines+1),
#      xlab = "Time", ylab = "sim$calcium level")
# j=1
# for(t in which(sim$amplitudes[idx[j+1],]>0)) { segments(t, j, t, j+0.2, col = sim$cluster_neurons[idx[nlines+j+1]] )}
# lines(sim$calcium[idx[nlines+j+1],]/(const)+j, col = sim$cluster_neurons[idx[nlines+j+1]])
# for(j in 2:(nlines)) {
#   for(t in which(sim$amplitudes[idx[nlines+j],]>0)) { segments(t, j, t, j+0.2, col = sim$cluster_neurons[idx[nlines+j]]) }
#   lines(sim$calcium[idx[nlines+j],]/(const)+j, col = sim$cluster_neurons[idx[nlines+j]])
# }
# par(mfrow=c(1,1))
#---------# #-----------# #---------# #---------# #-----------# #---------#
#---------# #-----------# #---------# #---------# #-----------# #---------#
