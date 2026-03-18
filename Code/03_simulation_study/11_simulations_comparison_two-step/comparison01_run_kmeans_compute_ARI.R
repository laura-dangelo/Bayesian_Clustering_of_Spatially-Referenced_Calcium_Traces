
library(coca)
library(mclust)
library(mcclust.ext)

replications_sim = 50


for(n_sim in 1:replications_sim){
  nameopen = paste0("10_simulations_sensitivity/simulated_data/simulated_data", n_sim,".RDS")
  sim = readRDS(nameopen)
  
  TT = sim$TT
  N = sim$n
  
  name = paste0("11_simulations_comparison_two-step/results/estimated_calcium_JW_simu", n_sim, ".RDS")
  calcium = readRDS(name)

  
  n_clus = 2:10
  cluster_kmeans = array(dim=c(N,N,length(n_clus)))
  estimated_clusters = matrix(0, length(n_clus),N)
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
    # heatmap(output_cl)
    estimated_clusters[k,] = minVI(output_cl)$cl
  }
  chooseK = maximiseSilhouette( cluster_kmeans, estimated_clusters, 10 )
  
  estimated_cluster_kmeans = estimated_clusters[chooseK$K,]
  name = paste0("11_simulations_comparison_two-step/results/estimates_clusterkmeans_simu", n_sim, ".RDS")
  saveRDS(estimated_cluster_kmeans, file = name)
  

}


ARIs_kmeans = numeric(replications_sim)

for(n_sim in 1:replications_sim){
  nameopen = paste0("10_simulations_sensitivity/simulated_data/simulated_data", n_sim,".RDS")
  sim = readRDS(nameopen)
  
  TT = sim$TT
  N = sim$n

  name = paste0("11_simulations_comparison_two-step/results/estimates_clusterkmeans_simu", n_sim, ".RDS")
  estimated_cluster_kmeans = readRDS(file = name)
  
  ARIs_kmeans[n_sim] = adjustedRandIndex(sim$cluster_neurons, estimated_cluster_kmeans)
  
}

name = paste0("11_simulations_comparison_two-step/output_RDS/ARIs_kmeans.RDS")
saveRDS(ARIs_kmeans, file = name)
