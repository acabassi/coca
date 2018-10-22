#' Consensus clustering with k-means
#'
#' This function allows to perform consensus clustering using the k-means clustering algorithm, for a fixed number
#' of clusters. We consider the number clusters K to be fixed, for simplicity.
#' @param data N X P data matrix
#' @param K number of clusters
#' @param B of iterations
#' @param pItem proportion of items sampled at each iteration
#' @param clMethod clustering method
#' @param distHC distance used for hierarchical clustering
#' @return The output is a consensus matrix, that is a symmetric matrix where the element in
#' position (i,j) corresponds to the proportion of times that items i and j have been clustered
#' together.
#' @author Alessandra Cabassi \email{ac2051@cam.ac.uk}
#' @references Monti, S., Tamayo, P., Mesirov, J. and Golub, T., 2003. Consensus clustering: a resampling-based method for
#' class discovery and visualization of gene expression microarray data. Machine learning, 52(1-2), pp.91-118.
#' @examples
#' ## Load one dataset with 300 observations, 2 variables, 6 clusters
#' data <- as.matrix(read.csv(system.file("extdata", "dataset1.csv", package = "coca"),
#' row.names = 1))
#'
#' ## Compute consensus clustering with K=6 clusters
#' cm <- consensusCluster(data, K = 6)
#' @export

consensusCluster = function(data, K, B = 100, pItem = 0.8, clMethod = "km",
                            distHC = "euclidean"){

  # Save number of observations
  N <- dim(data)[1]
  # Create vector of data indices that will be used to update the co-clustering matrix
  dataIndices <- c(1:N)
  # Initialise empty co-clustering matrix and an auxiliar
  coClusteringMatrix <- indicatorMatrix <- matrix(0, N, N)

  # For each step of the algorithm
  for(b in 1:B){
    # Sample a proportion pItem of the observations without replacement
    items <- sample(N, ceiling(N*pItem), replace = FALSE)

    uniqueData <- unique(data[items,])
    nUniqueDataPoints <- nrow(uniqueData)
    if(nUniqueDataPoints>K){
        # If the chosen clustering method is k-means
        if(clMethod == "km")
            # Apply k-means to the subsample and extract cluster labels
            cl <- stats::kmeans(data[items,], K, iter.max = 1000)$cluster

        # If the chosen clustering method is hierarchical clustering
        else if(clMethod == "hc"){
            # Calculate pairwise distances between observations
            distances <- stats::dist(data, method = distHC)
            # Apply hierarchical clustering to the subsample
            hClustering <- stats::hclust(distances, method = "average")
            clLabels <- stats::cutree(hClustering, K)
        }

        # Update matrix containing counts of number of times each pair has been sampled together
        indicatorMatrix <- indicatorMatrix + crossprod(t(as.numeric(dataIndices%in%items)))

        # For each cluster
        for(k in 1:K){
            # Update counts of number of times each pair has been put in the same cluster
            coClusteringMatrix[items,items] <- coClusteringMatrix[items,items] + crossprod(t(as.numeric(cl==k)))
        }
    }
  }

  if(!sum(indicatorMatrix)==0){
      # Normalise co-clustering matrix
      consensusMatrix <- coClusteringMatrix/indicatorMatrix
  }else{
      consensusMatrix <- indicatorMatrix
      warning(paste("Consensus matrix is empty for K =", K, "because there are less than", K,
              "distinct data points", sep = ""))
  }

  return(consensusMatrix)

}
