#' Expand matrix of clusters
#'
#' Expand matrix of cluster labels into matrix of clusters
#'
#' @param clLabels Matrix of cluster labels of size N x M
#' @param datasetNames Vector of cluster names of length M. Default is null
#' @return The output is a list containing the matrix of clusters `moc`, a vector containing
#' the dataset indicator `datasetIndicator` and a vector of expanded `datasetNames` of
#' length sum(K)
#' @author Alessandra Cabassi \email{ac2051@cam.ac.uk}
#' @export

expandMOC = function(clLabels, datasetNames = NULL){

  # Number of data points
  N <- dim(clLabels)[1]
  # Number of datasets
  M <- dim(clLabels)[2]

  # Number of clusters in each dataset
  K <- rep(NA, M)
  for(i in 1:M){
    K[i] <- length(unique(clLabels[,i]))
  }

  # If no names are provided for the datasets, assign names "1", "2", "3" and so on.
  if(is.null(datasetNames)){
    datasetNames <- as.character(1:M)
  }

  # Initialise objects containing the output
  moc <- matrix(NA, N, sum(K))
  rownames(moc) <- rownames(clLabels)
  datasetIndicator <- NULL
  datasetNamesMOC <- NULL

  count <- 0

  # For each dataset
  for(i in 1:M){
    # Get cluster names
    clLabels_i <- unique(clLabels[,i])
    # For each cluster name
    for(j in clLabels_i){
      # Fill in the corresponding row in MOC
      count <- count + 1
      moc[,count] <- (clLabels[i] == j)*1
      datasetIndicator <- c(datasetIndicator, i)
      datasetNamesMOC <- c(datasetNamesMOC, datasetNames[i])
    }
  }

  # Check whether all the rows of MOC have been filled
  if(count != dim(moc)[2])
    warning("moc has not been filled completely!")

  output <- list(moc = moc, datasetIndicator = datasetIndicator,
                 datasetNames = datasetNamesMOC)
  return(output)
}
