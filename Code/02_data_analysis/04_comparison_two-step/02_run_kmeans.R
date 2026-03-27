#---------# #-----------# #---------# #---------# #-----------# #---------# #---------# 
#---------#              STEP 2: CLUSTER THE DECONVOLVED SERIES             #---------#
#---------# #-----------# #---------# #---------# #-----------# #---------# #---------# 
#
# This script clusters the deconvolved series using the consensus kmeans algorithm contained in the package "coca"
#  _________________________________
#  YOU CAN AVOID RUNNING THIS SCRIPT 
#  ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
# the RDS files produced by this script are available in the "Code/02_data_analysis/04_comparison_two-step/output_RDS" folder.

# devtools::install_github("sarawade/mcclust.ext")
library(coca)
library(mclust)
library(mcclust.ext)

idx = readRDS("../Data/Time_windows/indices.RDS")

n_w = 0
n_window = idx[n_w+1]

idx = c(4,17,36,77)

for(n_window in idx) {
  n_w = n_w+1
  cat(paste("Running window",n_w,"out of",length(idx)," - ID:",n_window,"\n"))
  
  filename = paste0("../Data/Time_windows/calcium_window", n_window, ".RDS")
  calcium = readRDS(filename)
  
  TT = nrow(calcium)
  N = ncol(calcium)
  
  name = paste0("02_data_analysis/04_comparison_two-step/output_RDS/estimated_calcium_JW_window", n_window, ".RDS")
  estimated_calcium = readRDS(name)

  n_clus = 2:15
  cluster_kmeans = array(dim=c(N,N,length(n_clus)))
  estimated_clusters = matrix(0, length(n_clus),N)
  for(k in 1:length(n_clus)) {
    set.seed(1234)
    output_cl = consensusCluster(
      data = estimated_calcium,
      K = n_clus[k],
      B = 200,
      pItem = 0.8,
      clMethod = "kmeans",
      dist = "euclidean"
    )
    cluster_kmeans[,,k] = output_cl
    estimated_clusters[k,] = minVI(output_cl)$cl
  }
  chooseK = maximiseSilhouette( cluster_kmeans, estimated_clusters, maxK=15 )
  K = which.max(chooseK$silhouette[4:length(chooseK$silhouette)])+1
  
  estimated_cluster_kmeans = estimated_clusters[K,]
  name = paste0("02_data_analysis/04_comparison_two-step/output_RDS/estimates_clusterkmeans_window", n_window, ".RDS")
  saveRDS(estimated_cluster_kmeans, file = name)
  

}


