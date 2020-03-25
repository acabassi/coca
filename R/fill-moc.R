#' Fill Matrix-Of-Clusters
#'
#' This function fills in a matrix of clusters that contains NAs, by estimating
#' the missing cluster labels based on the available ones or based on the other
#' datasets. The predictive accuracy of this method can also be estimated via
#' cross-validation.
#'
#' @param clLabels N X M matrix containing cluster labels. Element (n,m)
#' contains the cluster label for element data point n in cluster m.
#' @param data List of M datasets to be used for the label imputation.
#' @param computeAccuracy Boolean. If TRUE, for each missing element, the
#' performance of the predictive model used to estimate the corresponding
#' missing label is computer. Default is FALSE.
#' @param verbose Boolean. If TRUE, for each NA, the size of the matrix used to
#' estimate its values is printed to screen. Default is FALSE.
#' @return The output is a list containing:
#' \item{fullClLabels}{the same matrix of clusters as the input matrix
#' \code{clLabels}, where NAs have been replaced by their estimates, where
#' possible.}
#' \item{nRows}{matrix where the item in position (i,j) indicates the
#' number of observations used in the predictive model used to estimate  the
#' corresponding missing label in the \code{fullClLabels} matrix.}
#' \item{nColumns}{matrix where the item in position (i,j) indicates the
#' number of covariates used in the predictive model used to
#' estimate the corresponding missing label in the \code{fullClLabels} matrix.}
#' \item{accuracy}{a matrix where each element
#' corresponds to the predictive accuracy of the predictive model used to
#' estimate the corresponding label in the cluster label matrix. This is only
#' returned if the argument \code{computeAccuracy} is set to TRUE.}
#' \item{accuracy_random}{This is computed in the same way as \code{accuracy},
#' but with the labels randomly shuffled. This can be used in order to assess
#' the predictive accuracy of the imputation algorithm and is returned only if
#' the argument \code{computeAccuracy} is set to TRUE.}
#' @author Alessandra Cabassi \email{alessandra.cabassi@mrc-bsu.cam.ac.uk}
#' @references The Cancer Genome Atlas, 2012. Comprehensive molecular portraits
#' of human breast tumours. Nature, 487(7407), pp.61–70.
#' @examples
#' # Load data
#' data <- list()
#' data[[1]] <- as.matrix(read.csv(system.file("extdata", "dataset1.csv",
#'                        package = "coca"), row.names = 1))
#' data[[2]] <- as.matrix(read.csv(system.file("extdata", "dataset2.csv",
#'                        package = "coca"), row.names = 1))
#' data[[3]] <- as.matrix(read.csv(system.file("extdata", "dataset3.csv",
#'                        package = "coca"), row.names = 1))
#'
#' # Build matrix of clusters
#' outputBuildMOC <- buildMOC(data, M = 3, K = 6, distances = "cor")
#'
#' # Extract matrix of clusters
#' clLabels <- outputBuildMOC$clLabels
#'
#' # Impute missing values using full datasets
#' outputFillMOC <- fillMOC(clLabels, data)
#'
#' # Extract full matrix of cluster labels
#' clLabels2 <- outputFillMOC$fullClLabels
#' @export

fillMOC <- function(clLabels,
                    data,
                    computeAccuracy = FALSE,
                    verbose = FALSE) {


    N <- dim(clLabels)[1]  # Number of data points
    M <- dim(clLabels)[2]  # Number of data sets

    # Save rownames
    allRowNames <- rownames(clLabels)

    # Remove data points that have NAs
    remove_rows <- allRowNames[unique(which(rowSums(is.na(clLabels)) > 0))]

    # Convert clLabels to data.frame
    clLabels <- as.data.frame(clLabels)
    row.names(clLabels) <- allRowNames

    # Initialise matrix of cluster labels with imputed values
    fullClLabels <- clLabels

    if (computeAccuracy) {
        # Initialise matrix containing accuracy for each element of the matrix
        # clLabels
        accuracy <- matrix(NA, N, M)
        accuracy_random <- matrix(NA, N, M)
    }

    # Initialise matrices containing number of rows and columns of each design
    # matrix used for prediction
    n_rows <- n_columns <- matrix(NA, N, M)

    if (is.null(data))
        stop("Please provide full datasets if you want them to be used for
                 the estimation.")

    # For each data point
    for (i in seq_len(N)) {
        label_i <- allRowNames[i]

        # For each dataset
        for (j in seq_len(M)) {
            # Check if label is missing
            if (is.na(clLabels[[j]][i])) {
                if (verbose)
                  print(paste("Element (", label_i, ",", names(clLabels)[j],
                              ") is NA", sep = ""))

                remove_rows_i <- remove_rows[!remove_rows == label_i]
                retain_rows_i <-
                    allRowNames[-which(allRowNames %in% remove_rows_i)]

                counter <- 0

                for (l in seq_len(M)) {
                  if (!is.na(clLabels[[l]][i])) {
                    counter <- counter + 1

                    if (counter == 1) {
                      # Design matrix for estimation
                      Xfit <- data[[l]][retain_rows_i, ]

                      # Desing matrix for prediction
                      Xpredict <- as.matrix(data[[l]][label_i, ])
                      Xpredict <- t(Xpredict)

                      # Response variable
                      response <- as.factor(clLabels[[j]][which(allRowNames %in%
                        retain_rows_i)])
                      names(response) <- retain_rows_i
                      response <- response[-which(names(response) == label_i)]
                      response <- as.factor(response)
                    } else {
                      # Append columns to design matrices
                      Xfit <- cbind(Xfit, data[[l]][retain_rows_i, ])
                      XpredictADD <- data[[l]][label_i, ]
                      XpredictADD <- t(XpredictADD)
                      Xpredict <- cbind(Xpredict, XpredictADD)
                    }
                  }
                }

                Xfit <- Xfit[-which(rownames(Xfit) == label_i), ]

                if (exists("Xfit")) {
                  # Fit glm
                  foldIDs <- stratifiedSamplingForCV(response)
                  glm_i <-
                      glmnet::cv.glmnet(Xfit, response, alpha = 1,
                                        foldid = foldIDs,
                                        family = "multinomial",
                                        type.multinomial = "grouped")

                  # Predict labels
                  prediction <- stats::predict(glm_i, newx = Xpredict,
                                               type = "class",
                    s = "lambda.min")
                  prediction <- as.vector(prediction)

                  if (verbose) {
                    print(paste("Predicted value for element (", i, " , ", j,
                                ") is: ", prediction, sep = ""))
                    print(paste("Matrix used to estimate element (", i, " , ",
                                j, ") is of size ", dim(Xfit)[1], " x ",
                                dim(Xfit)[2] - 1, sep = ""))
                  }

                  # Save number of rows and columns in desing matrix used to
                  # build the model
                  n_rows[i, j] <- dim(Xfit)[1]
                  n_columns[i, j] <- dim(Xfit)[2] - 1

                  ### Assess predictive performance through 5-fold CV
                  # Divide observations into 5 groups
                  if (computeAccuracy & length(table(response)) > 1) {
                    folds <- stratifiedSamplingForCV(response)
                    accuracy_l <- rep(NA, 5)
                    accuracy_random_l <- rep(NA, 5)
                    misclassRate_l <- rep(NA, 5)

                    for (l in seq_len(5)) {

                      XfitCV <- Xfit[-which(folds == l), ]
                      XpredCV <- Xfit[which(folds == l), ]

                      responseFitCV <- response[-which(folds == l)]
                      responsePredCV <- response[which(folds == l)]

                      foldIDs <- stratifiedSamplingForCV(responseFitCV)
                      glm_l <-
                          glmnet::cv.glmnet(XfitCV, responseFitCV,alpha = 1,
                                            foldid = foldIDs,
                                            family = "multinomial",
                                            type.multinomial = "grouped")
                      predictions_l <-
                          stats::predict(glm_l, newx = XpredCV, type = "class",
                                         s = "lambda.min")
                      predictions_l <- as.vector(predictions_l)

                      if (verbose) {
                        print(paste("CV predictions ", l, sep = ""))
                        print(predictions_l)
                      }

                      # Compute predictive accuracy
                      accuracy_l[l] <-
                          caret::postResample(predictions_l,
                                              responsePredCV)[[1]]

                      # Also generate random prediction
                      random_predictions_l <-
                          sample(responseFitCV, size = length(predictions_l),
                                 replace = FALSE)
                      # Compute predictive accuracy for random prediction
                      accuracy_random_l[l] <-
                          caret::postResample(random_predictions_l,
                                              responsePredCV)[[1]]
                    }

                    # Compute average accuracy
                    accuracy[i, j] <- mean(accuracy_l)
                    # Compute average accuracy of random prediction
                    accuracy_random[i, j] <- mean(accuracy_random_l)

                    if (verbose) {
                      print(paste("Prediction accuracy for element (", i, " , ",
                        j, ") is", accuracy[i, j], sep = ""))
                      print(paste("Random accuracy for element (", i, " , ", j,
                                  ") is", accuracy_random[i, j], sep = ""))
                    }
                  }

                  fullClLabels[i, j] <- prediction
                }
            }
        }
    }

    if (computeAccuracy) {
        output <-
            list(
                fullClLabels = fullClLabels,
                nRows = n_rows,
                nColumns = n_columns,
                accuracy = accuracy,
                accuracy_random = accuracy_random
            )
    } else {
        output <-
            list(fullClLabels = fullClLabels,
                 nRows = n_rows,
                 nColumns = n_columns)
    }

    return(output)
}

#' Divide data into 5 subsets using stratified sampling
#'
#' This function is used to do stratified subsampling based on the
#' number of observations in each group in the response
#'
#' @param response Vector of categorical responses
#' @return The function returns a vector of labels to assign each observation to
#' a different fold
#' @author Alessandra Cabassi  \email{alessandra.cabassi@mrc-bsu.cam.ac.uk}
#' @keywords internal
#'
stratifiedSamplingForCV <- function(response) {
    fold_labels <- rep(NA, length(response))

    for (g in unique(response)) {
        indices_g <- which(response == g)
        n_g <- length(indices_g)
        n_g_per_fold <- floor(n_g)/5

        for (h in seq_len(4)) {
            indices_g_fold_h <- sample(indices_g, n_g_per_fold)
            fold_labels[indices_g_fold_h] <- h
            indices_g <- indices_g[which(!(indices_g %in% indices_g_fold_h))]
        }

        fold_labels[indices_g] <- 5
    }

    return(fold_labels)
}
