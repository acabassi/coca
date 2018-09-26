#' Consensus clustering with k-means
#'
#' This function allows to perform consensus clustering using the k-means clustering algorithm, for a fixed number
#' of clusters. Contrarily to what was suggested in the original paper by Monti et al. (2003), here we consider the number
#' clusters K to be fixed, for simplicity.
#' @param data N X P data matrix
#' @param K number of clusters
#' @param B of iterations
#' @param pItem proportion of items sampled at each iteration
#' @param clMethod clustering method
#' @return The output is a consensus matrix, that is a symmetric matrix where the element in position (i,j) corresponds to
#' the proportion of times that items i and j have been clustered together.
#' @author Alessandra Cabassi \email{ac2051@cam.ac.uk}
#' @references Monti, S., Tamayo, P., Mesirov, J. and Golub, T., 2003. Consensus clustering: a resampling-based method for
#' class discovery and visualization of gene expression microarray data. Machine learning, 52(1-2), pp.91-118.
#' @export

consensusCluster = function(data, K, B = 100, pItem = 0.8, clMethod = "km", distHC = "euclidean"){

  N <- dim(data)[1]
  dataIndices <- c(1:N)
  coClusteringMatrix <- indicatorMatrix <- matrix(0, N, N)

  for(b in 1:B){
    items <- sample(N, ceiling(N*0.8), replace = FALSE)

    if(clMethod == "km")
      cl <- kmeans(data[items,], K, iter.max = 100)$cluster
    else if(clMethod == "hc"){
      distances <- dist(data, method = distHC)
      hClustering <- hclust(distances, method = "average")
      clusterLabels <- cutree(hClustering, K)
    }

    indicatorMatrix <- indicatorMatrix + crossprod(t(as.numeric(dataIndices%in%items)))
    for(k in 1:K){
      coClusteringMatrix[items,items] <- coClusteringMatrix[items,items] + crossprod(t(as.numeric(cl==k)))
    }
  }
  coClusteringMatrix/indicatorMatrix
}
