#' Build Matrix-Of-Clusters
#'
#' This function creates a matrix of clusters, starting from a list of heterogeneous datasets.
#'
#' @param data List of M datasets, each of size N X P_m, m = 1, ..., M.
#' @param M number of datasets.
#' @param K Vector containing the number of clusters in each dataset. If given an integer
#' instead of a vector it is assumed that each dataset has the same number of clusters
#' @param methods A vector of strings containing the names of the clustering methods to be
#' used to cluster the observations in each dataset
#' @param full Boolean. If TRUE, if there are any missing observations in one or more datasets,
#' the corresponding cluster labels will be estimated through generalised linear models on the
#' basis of the available labels.
#' @return The output is the Matrix-Of-Clusters (MOC). It is a binary matrix of size N x sum(K)
#' where element (n,k) contains a 1 if observation n belongs to the corresponding cluster,
#' 0 otherwise. It also returns a vector called datasetIndicator of length sum(K) in which
#' each element is the number of the dataset to which the cluster belongs.
#' @author Alessandra Cabassi \email{ac2051@cam.ac.uk}
#' @references The Cancer Genome Atlas, 2012. Comprehensive molecular portraits of human
#' breast tumours. Nature, 487(7407), pp.61–70.
#' @export

buildMOC = function(data, M, K, methods, full = FALSE){

  ### Match data ###
  obsNames <- rownames(data[[1]])
  for(i in 2:M){
      obsNames <- unique(c(obsNames, rownames(data[[i]])))
  }
  N <- length(obsNames)

  ### Sum number of clusters for each dataset ###
  if(length(K)==M){
    Ktot = sum(K)
  }else if(length(K)==1){
    Ktot = K*M
    K = rep(K, M)
  }else{
    stop("K must be either a vector of length M or a scalar.")
  }

  # Initalise empty matrices
  moc = matrix(NA, N, Ktot) # Binary matrix
  clLabels = matrix(NA, N, M) # Matrix containing cluster numbers
  datasetIndicator <- rep(NA, Ktot)
  # (it is equivalent to the matrix of clusters, but will be used to fill in NAs more quickly)

  rownames(moc) <- obsNames
  rownames(clLabels) <- obsNames

  count_k <- 0

  if(length(methods)==1)
    methods = rep(methods, M)
  # Should add another check here, on the length of vector methods

  # For each dataset
  for(i in 1:M){

    # Choose clustering algorithm
    method_i <- methods[i]

    # If it is k-means
    if(method_i == "kmeans"){

      # Find cluster labels
      newClusterLabels <-  kmeans(data[[i]], K[i])$cluster
      clLabels[names(newClusterLabels),i] <- newClusterLabels

      # Store them in moc matrix
      for(j in unique(newClusterLabels)){
        count_k <- count_k + 1
        bin <- rep(NA, N)
        names(bin) <- rownames(moc)
        bin[names(newClusterLabels)] <- (newClusterLabels == j)*1
        moc[,count_k] = bin
        datasetIndicator[count_k] <- i
      }

     }
    else if(method_i == "hclust"){

      # Compute distances between data points in dataset i
      d <- stats::as.dist(1 - cor(t(data[[i]]), method = "pearson"))

      # Find clusters through hierarchical clustering
      hCl <- stats::hclust(d, method = "average")

      # Extract cluster labels
      newClusterLabels <-  stats::cutree(hCl, K[i])
      clLabels[names(newClusterLabels),i] <- newClusterLabels

      # Store them in moc matrix
      for(j in unique(newClusterLabels)){
        count_k <- count_k + 1
        bin <- rep(NA, N)
        names(bin) <- rownames(moc)
        bin[names(newClusterLabels)] <- (newClusterLabels == j)*1
        moc[,count_k] = bin
        datasetIndicator[count_k] <- i
      }

    }else{
      stop("Clustering method name not recognised.")
    }
  }

  if(count_k != Ktot){
    stop("Something went wrong: matrix of clusters has not been filled properly.")
  }

  number_nas = sum(is.na(moc))

  # Should return also dataset indicator
  output <- list(moc = moc, datasetIndicator = datasetIndicator,
                 number_nas = number_nas, clLabels = clLabels)
  return(output)
}
