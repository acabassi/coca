#' Fill Matrix-Of-Clusters
#'
#' This function fills in a matrix of clusters that contains NAs, by estimating the missing cluster labels
#' based on the available ones.
#'
#' @param clLabels N X M matrix containing cluster labels. Element (n,m) contains the cluster label for element
#' @return The output is the same matrix of clusters, where NAs have been replaced by their estimates, where possible.
#' @author Alessandra Cabassi \email{ac2051@cam.ac.uk}
#' @references The Cancer Genome Atlas, 2012. Comprehensive molecular portraits of human breast tumours. Nature,
#' 487(7407), pp.61–70.
#' @export

fillMOC <- function(clLabels){

  N <- dim(clLabels)[1] # Number of datapoints
  M <- dim(clLabels)[2] # Number of datasets

  # clLabels <- matrix(as.character(clLabels), nrow = N, ncol = M, byrow = TRUE)
  clLabels <- as.data.frame(clLabels)

  fullClLabels <- clLabels

  # For each dataset
  for(i in 1:M){
    print(paste("M = ", i, sep = ""))

    # Check which observations are missing in dataset i
    nas_dataset_i <- is.na(clLabels[[i]])

    # If there is at least one observation missing in dataset i
    if(sum(nas_dataset_i)>0){
      print(paste("Column ", i, " has ", sum(nas_dataset_i)," NA(s)", sep = ""))

      # Build design matrix
      X <- data.frame(response = as.factor(clLabels[[i]]))
      count <- 1
      for(j in 1:M){
        # Only include columns corresponding to datasets for which there are no NAs
        # for the observations for which we want to estimate the i-th cluster label
        if(sum(is.na(clLabels[[j]][nas_dataset_i]))==0){
          X <- cbind(X, clLabels[[j]])
          count <- count + 1
          names(X)[count] <- names(clLabels)[j]
        }
      }

      if(!is.null(X)){

        # Fit glm
        Xfit <- X[which(!nas_dataset_i),]
        glm_i <- nnet::multinom(response ~ ., data=Xfit)

        # Predict labels
        Xpredict <- X[which(nas_dataset_i),]
        predictions <- stats::predict(glm_i, newdata = Xpredict, type = "class")
        print(paste("predicted values for NAs in column ", i, " are:", sep = ""))
        print(predictions)

        fullClLabels[nas_dataset_i,i] <- predictions
      } else {
        warning(paste("It was not possible to estimate all cluster labels for dataset",
                      i, sep = ""))
      }
    }else{
      print(paste("Column ", i, " does not have any NAs", sep = ""))
    }

  }
  return(fullClLabels)
}
