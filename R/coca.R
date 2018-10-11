#' Cluster-Of-Clusters Analysis
#'
#' This function allows to do Cluster-Of-Clusters-Analysis on a binary matrix where each column is a clustering of the data,
#' each row corresponds to a data point and the element in position (i,j) is equal to 1 if data point i belongs to cluster
#' j, 0 otherwise.
#'
#' @param MOC N X C data matrix, where C is the total number of clusters considered.
#' @param K Number of clusters.
#' @param B Number of iterations of the Consensus Clustering step.
#' @param pItem Proportion of items sampled at each iteration of the Consensus Cluster step.
#' @param hclustMethod Agglomeration method to be used by the hclust function for the hierarchical clustering step. Can be "single", "complete", "average", etc. For more details please see ?hclust.
#' @return The output is a consensus matrix, that is a symmetric matrix where the element in position (i,j) corresponds to
#' the proportion of times that items i and j have been clustered together and a vector of cluster labels.
#' @author Alessandra Cabassi \email{ac2051@cam.ac.uk}
#' @references The Cancer Genome Atlas, 2012. Comprehensive molecular portraits of human breast tumours. Nature,
#' 487(7407), pp.61–70.
#' @examples
#' ## Load data
#' data <- list()
#' data[[1]] <- as.matrix(read.csv(system.file("extdata", "dataset1.csv",
#' package = "coca"), row.names = 1))
#' data[[2]] <- as.matrix(read.csv(system.file("extdata", "dataset2.csv",
#' package = "coca"), row.names = 1))
#' data[[3]] <- as.matrix(read.csv(system.file("extdata", "dataset3.csv",
#' package = "coca"), row.names = 1))
#'
#' ## Set constants
#' n_datasets <- 3
#' n_clusters <- 6
#' N <- 300
#' true_labels <- as.matrix(read.csv(system.file("extdata", "cluster_labels.csv",
#' package = "coca"), row.names = 1))
#'
#' ## Fill label matrix with clusterings found with the k-means clustering algorithm
#' labelMatrix <- array(NA, c(n_clusters, N, n_datasets))
#' for(i in 1:n_datasets){
#'   output <- kmeans(data[[i]], n_clusters)
#'   for(k in 1:n_clusters){
#'     labelMatrix[k,,i] <- (output$cluster==k)
#'   }
#' }
#' # Convert label matrix from logic to numeric matrix
#' labelMatrix <- labelMatrix*1
#' @export

coca = function(MOC, K, B = 1000, pItem = 0.8, hclustMethod = 'average'){

  # Intialise output list
  output = list()

  ### Step 1. Compute the consensus matrix ###
  output$consensusMatrix <- consensusCluster(MOC, K, B, pItem)

  ### Step 2. Use hierarchical clustering on the consensus matrix ###
  distances <- stats::as.dist(1 - output$consensusMatrix)
  hClustering <- stats::hclust(distances, method = hclustMethod)
  output$clusterLabels <- stats::cutree(hClustering, K)

  return(output)

}
