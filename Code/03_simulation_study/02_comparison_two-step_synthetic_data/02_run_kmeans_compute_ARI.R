#-------------------# #-------------------# #-------------------#
#              CLUSTER SERIES USING CONSENSUS K-MEANS
#-------------------# #-------------------# #-------------------#

# These scripts replicate the comparison with the two-step approach reported in Section S4.3 of the Supplementary Material.
# Specifically, this script uses the consensus k-means algorithm of the package coca to obtain the clusters.

library(coca)
library(mclust)
library(mcclust.ext)

replications_sim = 50


for(n_sim in 1:replications_sim){
  nameopen = paste0("03_simulation_study/01_sensitivity_study/simulated_data/simulated_data", n_sim,".RDS")
  sim = readRDS(nameopen)
  
  TT = sim$TT
  N = sim$n
  
  name = paste0("03_simulation_study/02_comparison_two-step_synthetic_data/results/estimated_calcium_JW_simu", n_sim, ".RDS")
  calcium = readRDS(name)

  
  n_clus = 2:10
  cluster_kmeans = array(dim=c(N,N,length(n_clus)))
  estimated_clusters = matrix(0, length(n_clus),N)
  set.seed(123)
  for(k in 1:length(n_clus)) {
    output_cl = consensusCluster(
      data = calcium,
      K = n_clus[k],
      B = 200,
      pItem = 0.5,
      clMethod = "kmeans",
      dist = "euclidean"
    )
    cluster_kmeans[,,k] = output_cl
    estimated_clusters[k,] = minVI(output_cl)$cl
  }
  chooseK = maximiseSilhouette( cluster_kmeans, estimated_clusters, 10 )
  
  estimated_cluster_kmeans = estimated_clusters[chooseK$K,]
  name = paste0("03_simulation_study/02_comparison_two-step_synthetic_data/results/estimates_clusterkmeans_simu", n_sim, ".RDS")
  saveRDS(estimated_cluster_kmeans, file = name)
  

}


ARIs_kmeans = numeric(replications_sim)

for(n_sim in 1:replications_sim){
  nameopen = paste0("03_simulation_study/01_sensitivity_study/simulated_data/simulated_data", n_sim,".RDS")
  sim = readRDS(nameopen)
  
  TT = sim$TT
  N = sim$n

  name = paste0("03_simulation_study/02_comparison_two-step_synthetic_data/results/estimates_clusterkmeans_simu", n_sim, ".RDS")
  estimated_cluster_kmeans = readRDS(file = name)
  
  ARIs_kmeans[n_sim] = adjustedRandIndex(sim$cluster_neurons, estimated_cluster_kmeans)
  
}

name = paste0("03_simulation_study/02_comparison_two-step_synthetic_data/output_RDS/ARIs_kmeans.RDS")
saveRDS(ARIs_kmeans, file = name)
