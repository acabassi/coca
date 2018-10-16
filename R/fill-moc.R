#' Fill Matrix-Of-Clusters
#'
#' This function fills in a matrix of clusters that contains NAs, by estimating the missing cluster labels
#' based on the available ones.
#'
#' @param clLabels N X M matrix containing cluster labels. Element (n,m) contains the cluster label
#' for element data point n in cluster m.
#' @param computeAccuracy Boolean. If TRUE, for each missing element, the performance of the
#' predictive model used to estimate the corresponding missing label is computer.
#' @param verbose Boolean. If TRUE, for each NA, the size of the matrix used to estimate its values is printed
#' to screen
#' @return The output is the same matrix of clusters, where NAs have been replaced by their estimates,
#' where possible. If `computeAccuracy` is TRUE, then also an object called `accuracy` is returned,
#' where each
#' @author Alessandra Cabassi \email{ac2051@cam.ac.uk}
#' @references The Cancer Genome Atlas, 2012. Comprehensive molecular portraits of human breast tumours.
#' Nature, 487(7407), pp.61–70.
#' @export

fillMOC <- function(clLabels, computeAccuracy = FALSE, verbose = FALSE){

  N <- dim(clLabels)[1] # Number of datapoints
  M <- dim(clLabels)[2] # Number of datasets

  # clLabels <- matrix(as.character(clLabels), nrow = N, ncol = M, byrow = TRUE)
  clLabels <- as.data.frame(clLabels)

  # Initialise matrix of cluster labels with imputed values
  fullClLabels <- clLabels

  if(computeAccuracy){
    # Initialise matrix containing accuracy for each element of the matrix clLabels
    accuracy <- matrix(NA, N, M)
  }
  n_rows <- n_columns <- matrix(NA, N, M)

  # For each data point
  for(i in 1:N){

    # For each dataset
    for(j in 1:M){

      # Check if label is missing
      if(is.na(clLabels[[j]][i])){

        if(verbose)
          print(paste("Element", i, j, "is NA", sep =" "))

          # Build design matrix
          X <- data.frame(response = as.factor(clLabels[[j]]))
          count <- 1
            for(l in 1:M){
              if(!is.na(clLabels[[l]][i])){
                X <- cbind(X, clLabels[[l]])
                count <- count + 1
                names(X)[count] <- names(clLabels)[l]
              }
            }

          if(!is.null(X)){

            # Fit glm
            remove_rows <- unique(c(which(rowSums(is.na(clLabels))>0),i))

            Xfit <- X[-remove_rows,]
            n_rows[i,j] <- dim(Xfit)[1]
            n_columns[i,j] <- dim(Xfit)[2]-1

            glm_i <- nnet::multinom(response ~ ., data=Xfit)

            # Predict labels
            Xpredict <- X[i,]
            prediction <- stats::predict(glm_i, newdata = Xpredict, type = "class")


            if(verbose) {
              print(paste("Predicted value for element (", i, " , ", j, ") is: ", prediction, sep = ""))
              print(paste("Matrix used to estimate element (", i, " , ", j, ") is of size ", dim(Xfit)[1],
                          " x ", dim(Xfit)[2]-1, sep = ""))
            }

            ### Assess predictive performance through 5-fold CV
            # Divide observations into 5 groups
            if(computeAccuracy & length(table(X$response))>1){
              folds <- stratifiedSamplingForCV(Xfit$response)
              accuracy_l <- rep(NA, 5)
              for(l in 1:5){
                XfitCV <- Xfit[-which(folds==l),]
                glm_l <- nnet::multinom(response ~ ., data = XfitCV)
                XpredCV <- Xfit[which(folds==l),]
                predictions_l <- stats::predict(glm_l, newdata = XpredCV, type ="class")
                accuracy_l[l] <- caret::postResample(predictions_l, XpredCV$response)[[1]]
              }
              accuracy[i,j] <- mean(accuracy_l)
            }
            if(verbose)
              print(paste("Prediction accuracy for element (", i, " , ", j, ") is", accuracy[i,j], sep = ""))

            fullClLabels[i,j] <- prediction
        }else{
          warning(paste("It was not possible to estimate the cluster label for datapoint", i,"
                        dataset ",j, sep = ""))
        }
      }
    }
  }

  if(computeAccuracy){
    output <- list(fullClLabels = fullClLabels, nRows = n_rows, nColumns = n_columns, accuracy = accuracy)
  }else{
    output <- list(fullClLabels = fullClLabels, nRows = n_rows, nColumns = n_columns)
  }

  return(output)
}

#' Divide data into 5 subsets using stratified sampling
#'
#' This function is used to do stratified subsampling based on the
#' number of observations in each group in the response
#'
#' @param response Vector of categorical responses
#' @return The function returns a vector of labels to assign each observation to a different fold
#' @author Alessandra Cabassi  \email{ac2051@cam.ac.uk}
#'
stratifiedSamplingForCV <- function(response){
  fold_labels <- rep(NA, length(response))
  for(g in unique(response)){
    indices_g <- which(response==g)
    n_g <- length(indices_g)
    n_g_per_fold <- floor(n_g)/5
    for(h in 1:4){
      indices_g_fold_h <- sample(indices_g, n_g_per_fold)
      fold_labels[indices_g_fold_h] <- h
      indices_g <- indices_g[which(!(indices_g %in% indices_g_fold_h))]
    }
    fold_labels[indices_g] <- 5
  }
  return(fold_labels)
}
