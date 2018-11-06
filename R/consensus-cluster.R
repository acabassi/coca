#' Consensus clustering with k-means
#'
#' This function allows to perform consensus clustering using the k-means clustering algorithm,
#' for a fixed number of clusters. We consider the number of clusters K to be fixed.
#' @param data N X P data matrix
#' @param K Number of clusters.
#' @param B Number of iterations.
#' @param pItem Proportion of items sampled at each iteration.
#' @param clMethod Clustering algorithm. Can be "hc" for hierarchical clustering, "km"
#' for k-means clustering, or "pam" for partitioning around medoids. Default is "hc". However,
#' if the data contain at least one covariate that is a factor, the default clustering
#' algorithm is "pam".
#' @param distHC Distance used for hierarchical clustering. Can be "pearson" (for 1 - Pearson
#' correlation), "spearman" (for 1- Spearman correlation), or any of the distances provided in
#' stats::dist() (i.e. "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski").
#' @return The output is a consensus matrix, that is a symmetric matrix where the element in
#' position (i,j) corresponds to the proportion of times that items i and j have been clustered
#' together.
#' @author Alessandra Cabassi \email{ac2051@cam.ac.uk}
#' @references Monti, S., Tamayo, P., Mesirov, J. and Golub, T., 2003. Consensus clustering:
#' a resampling-based method for class discovery and visualization of gene expression microarray
#' data. Machine learning, 52(1-2), pp.91-118.
#' @examples
#' ## Load one dataset with 300 observations, 2 variables, 6 clusters
#' data <- as.matrix(read.csv(system.file("extdata", "dataset1.csv", package = "coca"),
#' row.names = 1))
#'
#' ## Compute consensus clustering with K=6 clusters
#' cm <- consensusCluster(data, K = 6)
#' @export

consensusCluster = function(data, K, B = 100, pItem = 0.8, clMethod = "hc",
                            distHC = "euclidean"){

  # Save number of observations
  N <- dim(data)[1]
  # Save number of covariates
  P <- dim(data)[2]

  # Create vector of data indices that will be used to update the co-clustering matrix
  dataIndices <- c(1:N)
  # Initialise empty co-clustering matrix and an auxiliar
  coClusteringMatrix <- indicatorMatrix <- matrix(0, N, N)

  # Check wheter there are factors among the covariates
  containsFactors <- 0
  for(i in 1:P){
    containsFactors <- as.numeric(is.factor(data[,i]))+containsFactors
  }

  # For each step of the algorithm
  for(b in 1:B){
    # Sample a proportion pItem of the observations without replacement
    items <- sample(N, ceiling(N*pItem), replace = FALSE)

    # If there are more unique data points than clusters
    uniqueData <- unique(data[items,])
    nUniqueDataPoints <- nrow(uniqueData)
    if(nUniqueDataPoints>K){
        # If the chosen clustering method is PAM or there is at least one
        # covariate that is a factor
        if(clMethod == "pam" | containsFactors){
            # Apply pam to the subsample and extract cluster labels
            cl <- cluster::pam(data[items,], K)$clustering
        }
        # If the chosen clustering method is k-means
        else if(clMethod == "km"){
            # Apply k-means to the subsample and extract cluster labels
            cl <- stats::kmeans(data[items,], K, iter.max = 1000)$cluster
        }
        # If the chosen clustering method is hierarchical clustering
        else if(clMethod == "hc"){
            if(distHC == 'pearson' | distHC == 'spearman'){
                pearsonCor <- stats::cor(t(data[items,]), method = distHC)
                distances <- stats::as.dist(1-pearsonCor)
            }else{
                # Calculate pairwise distances between observations
                distances <- stats::dist(data[items,], method = distHC)
            }

            # Apply hierarchical clustering to the subsample
            hClustering <- stats::hclust(distances, method = "average")
            cl <- stats::cutree(hClustering, K)
        }else{
            stop('Clustering algorithm name not recognised. Please choose either `km` (for k-means
  clustering or `hc` for hierarchical clustering.')
        }

        # Update matrix containing counts of number of times each pair has been sampled together
        indicatorMatrix <- indicatorMatrix + crossprod(t(as.numeric(dataIndices%in%items)))

        # For each cluster
        for(k in 1:K){
            # Update counts of number of times each pair has been put in the same cluster
            coClusteringMatrix[items,items] <- coClusteringMatrix[items,items] +
                crossprod(t(as.numeric(cl==k)))
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
