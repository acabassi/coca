#' Build Matrix-Of-Clusters
#'
#' This function creates a matrix of clusters starting from a list of
#' heterogeneous datasets.
#'
#' @param data List of M datasets, each of size N X P_m, where m = 1, ..., M.
#' @param M Number of datasets.
#' @param K Vector containing the number of clusters in each dataset. If given
#' an integer instead of a vector it is assumed that each dataset has the same
#' number of clusters. If NULL, it is assumed that the true cluster numbers are
#' not known, therefore they will be estimated using the silhouette method.
#' @param maxK Vector of maximum cluster numbers to be considered for each
#' dataset if K is NULL. If given an integer instead of a vector it is assumed
#' that for each dataset the same maximum number of clusters must be considered.
#' Default is 10.
#' @param methods Vector of strings containing the names of the clustering
#' methods to be used to cluster the observations in each dataset. Each can be
#' "kmeans" (k-means clustering), "hclust" (hierarchical clustering), or "pam"
#' (partitioning around medoids). If the vector is of length one, the same
#' clustering method is applied to all the datasets. Default is "hclust".
#' @param distances Distances to be used in the clustering step for each
#' dataset. If only one string is provided, then the same distance is used for
#' all datasets. If the number of strings provided is the same as the number of
#' datasets, then each distance will be used for the corresponding dataset.
#' Default is "euclidean". Please note that not all distances are compatible
#' with all clustering methods. "euclidean" and "manhattan" work with all
#' available clustering algorithms. "gower" distance is only available for
#' partitioning around medoids. In addition, "maximum", "canberra", "binary" or
#' "minkowski" are available for k-means and hierarchical clustering.
#' @param fill Boolean. If TRUE, if there are any missing observations in one or
#' more datasets, the corresponding cluster labels will be estimated through
#' generalised linear models on the basis of the available labels.
#' @param computeAccuracy Boolean. If TRUE, for each missing element, the
#' performance of the predictive model used to estimate the corresponding
#' missing label is computer.
#' @param fullData Boolean. If TRUE, the full data matrices are used to estimate
#' the missing cluster labels (instead of just using the cluster labels of the
#' corresponding datasets).
#' @param savePNG Boolean. If TRUE, plots of the silhouette for each datasets
#' are saved as png files. Default is FALSE.
#' @param fileName If \code{savePNG} is TRUE, this is the string containing the
#' name of the output files. Can be used to specify the folder path too. Default
#' is "buildMOC". The ".png" extension is automatically added to this string.
#' @param widestGap Boolean. If TRUE, compute also widest gap index to choose
#' best number of clusters. Default is FALSE.
#' @param dunns Boolean. If TRUE, compute also Dunn's index to choose best
#' number of clusters. Default is FALSE.
#' @param dunn2s Boolean. If TRUE, compute also alternative Dunn's index to
#' choose best number of clusters. Default is FALSE.
#' @return This function returns a list containing:
#' \item{moc}{the Matrix-Of-Clusters, a binary matrix of size N x sum(K)
#' where element (n,k) contains a 1 if observation n belongs to the
#' corresponding cluster, 0 otherwise.}
#' \item{datasetIndicator}{a vector of length sum(K) in which
#' each element is the number of the dataset to which the cluster belongs.}
#' \item{number_nas}{the total number of NAs in the matrix of clusters. (If the
#' MOC has been filled with imputed values, \code{number_nas} indicates the
#' number of NAs in the original MOC.)}
#' \item{clLabels}{a matrix that is equivalent to the matrix of clusters, but is
#' in compact form, i.e. each column corresponds to a dataset, each row
#' represents an observation, and its values indicate the cluster labels.}
#' \item{K}{vector of cluster numbers in each dataset. If these are provided as
#' input, this is the same as the input (expanded to a vector if the input is an
#' integer). If the cluster numbers are not provided as input, this vector
#' contains the cluster numbers chosen via silhouette for each dataset.}
#' @author Alessandra Cabassi \email{alessandra.cabassi@mrc-bsu.cam.ac.uk}
#' @references The Cancer Genome Atlas, 2012. Comprehensive molecular portraits
#' of human breast tumours. Nature, 487(7407), pp.61–70.
#' @references Rousseeuw, P.J., 1987. Silhouettes: a graphical aid to the
#' interpretation and validation of cluster analysis. Journal of computational
#' and applied mathematics, 20, pp.53-65.
#' @examples
#' # Load data
#' data <- list()
#' data[[1]] <- as.matrix(read.csv(system.file("extdata", "dataset1.csv",
#' package = "coca"), row.names = 1))
#' data[[2]] <- as.matrix(read.csv(system.file("extdata", "dataset2.csv",
#' package = "coca"), row.names = 1))
#' data[[3]] <- as.matrix(read.csv(system.file("extdata", "dataset3.csv",
#' package = "coca"), row.names = 1))
#'
#' # Build matrix of clusters
#' outputBuildMOC <- buildMOC(data, M = 3, K = 6, distances = "cor")
#'
#' # Extract matrix of clusters
#' matrixOfClusters <- outputBuildMOC$moc
#'
#' @export

buildMOC <-
    function(data,
             M,
             K = NULL,
             maxK = 10,
             methods = "hclust",
             distances = "euclidean",
             fill = FALSE,
             computeAccuracy = FALSE,
             fullData = FALSE,
             savePNG = FALSE,
             fileName = "buildMOC",
             widestGap = FALSE,
             dunns = FALSE,
             dunn2s = FALSE) {

    ###### Match data ######
    obsNames <- rownames(data[[1]])
    if (is.null(obsNames)) {
        N <- dim(data[[1]])[1]
        for (i in seq_len(M)) {
            if (is.null(dim(data[[i]])[1]) | !(dim(data[[i]])[1]) == N)
                stop("If datasets are not provided with rownames, then they need
                to have all the same number of rows (and the elements must be in
                the same order in each dataset), otherwise it is impossible
                to match observations from different dataset into the same
                matrix of clusters.")
            obsNames <- paste("Obs.", seq_len(N), sep = " ")
            rownames(data[[i]]) <- obsNames
        }
    } else {
        for (i in 2:M) {
            obsNames <- unique(c(obsNames, rownames(data[[i]])))
        }
        N <- length(obsNames)
    }

    if (length(methods) == 1) {
        methods <- rep(methods, M)
    } else if (length(methods) != M) {
        stop("The lenght of vector 'methods' must be either 1 (same clustering
            algorithm applied to all datasets) or M.")
    }

    if (length(distances) == 1) {
        distances <- rep(distances, M)
    } else if (length(distances) != M) {
        stop("The lenght of vector 'distances' must be either 1 (same clustering
        algorithm applied to all datasets) or M.")
    }

    ###### Choose cluster numbers ######
    if (is.null(K)) {
        if (length(maxK) == 1) {
            maxK <- rep(maxK, M)
        }

        K <- rep(NA, M)

        # For each dataset
        for (i in seq_len(M)) {
            N_i <- dim(data[[i]])[1]
            maxK_i <- maxK[i]
            tempClLabels <- matrix(NA, maxK_i - 1, N_i)
            distanceMatrix <- array(NA, c(N_i, N_i, maxK_i))

            # Choose clustering algorithm
            method_i <- methods[i]
            distance_i <- distances[i]

            # If clustering algorithm is kmeans
            if (method_i == "kmeans") {
                for (j in 2:maxK_i) {
                    ### Step 1. Compute the distance matrix
                    distanceMatrix[, , j - 1] <- as.matrix(
                        stats::dist(data[[i]], method = distance_i))
                    ### Step 2. Use k-means clustering to find cluster labels
                    tempClLabels[j - 1, ] <- stats::kmeans(data[[i]], j)$cluster
                }
            } else if (method_i == "hclust") {
                for (j in 2:maxK_i) {
                  ### Step 1. Compute the distance matrix
                  distanceMatrix[, , j - 1] <- as.matrix(
                      stats::dist(data[[i]], method = distance_i))
                  ### Step 2. Use hierarchical clust. on the consensus matrix
                  hClustering <- stats::hclust(
                      stats::dist(data[[i]], method = distance_i),
                      method = "complete")
                  tempClLabels[j - 1, ] <- stats::cutree(hClustering, j)
                }

                # If clustering algorithm is none of the above, stop.
            } else if (method_i == "pam") {
                for (j in 2:maxK_i) {
                  ### Step 1. Compute the distance matrix
                  if (distance_i == "cor") {
                    distanceMatrix[, , j - 1] <-
                        as.matrix(1 - stats::cor(t(data[[i]])))
                  } else {
                    distanceMatrix[, , j - 1] <-
                        as.matrix(cluster::daisy(data[[i]],
                      metric = distance_i))
                  }
                  ### Step 2. Use k-means clustering to find cluster labels
                  tempClLabels[j - 1, ] <-
                      cluster::pam(stats::as.dist(
                          distanceMatrix[,, j - 1]), j)$clustering
                }
            } else {
                stop("Clustering method name not recognised.")
            }
            K[i] <- maximiseSilhouette(
                distanceMatrix, tempClLabels, maxK_i, savePNG = savePNG,
                fileName = paste(fileName, "_dataset", i, sep = ""),
                isDistance = TRUE, widestGap = widestGap, dunns = dunns,
                dunn2s = dunn2s)$K
        }
    }

    ###### Sum number of clusters for each dataset ######
    if (length(K) == M) {
        Ktot <- sum(K)
    } else if (length(K) == 1) {
        Ktot <- K * M
        K <- rep(K, M)
    } else {
        stop("K must be either a vector of length M or a scalar.")
    }

    ###### Produce matrix of clusters #####
    # Initalise empty matrices
    moc <- matrix(NA, N, Ktot)  # Binary matrix
    clLabels <- matrix(NA, N, M)  # Matrix containing cluster numbers
    datasetIndicator <- rep(NA, Ktot)
    # (it is equivalent to the matrix of clusters, but will be used to fill in
    # NAs more quickly)

    rownames(moc) <- obsNames
    rownames(clLabels) <- obsNames

    count_k <- 0

    # For each dataset
    for (i in seq_len(M)) {
        # Choose clustering algorithm
        method_i <- methods[i]

        # If it is k-means
        if (method_i == "kmeans") {
            # Find cluster labels

            newClusterLabels <- stats::kmeans(data[[i]], K[i])$cluster
            clLabels[names(newClusterLabels), i] <- newClusterLabels

            # Store them in moc matrix
            for (j in unique(newClusterLabels)) {
                count_k <- count_k + 1
                bin <- rep(NA, N)
                names(bin) <- rownames(moc)
                bin[names(newClusterLabels)] <- (newClusterLabels == j) * 1
                moc[, count_k] <- bin
                datasetIndicator[count_k] <- i
            }
        } else if (method_i == "hclust") {
            # Compute distances between data points in dataset i
            d <- stats::dist(data[[i]], method = "euclidean")

            # Find clusters through hierarchical clustering
            hCl <- stats::hclust(d, method = "complete")

            # Extract cluster labels
            newClusterLabels <- stats::cutree(hCl, K[i])
            clLabels[names(newClusterLabels), i] <- newClusterLabels

            # Store them in moc matrix
            for (j in unique(newClusterLabels)) {
                count_k <- count_k + 1
                bin <- rep(NA, N)
                names(bin) <- rownames(moc)
                bin[names(newClusterLabels)] <- (newClusterLabels == j) * 1
                moc[, count_k] <- bin
                datasetIndicator[count_k] <- i
            }

            # If clustering algorithm is none of the above, stop.
        } else if (method_i == "pam") {
            # Find cluster labels

            newClusterLabels <- cluster::pam(data[[i]], K[i])$clustering
            clLabels[names(newClusterLabels), i] <- newClusterLabels

            # Store them in moc matrix
            for (j in unique(newClusterLabels)) {
                count_k <- count_k + 1
                bin <- rep(NA, N)
                names(bin) <- rownames(moc)
                bin[names(newClusterLabels)] <- (newClusterLabels == j) * 1
                moc[, count_k] <- bin
                datasetIndicator[count_k] <- i
            }

        } else {
            stop("Clustering method name must be either kmeans or hclust.")
        }
    }

    if (count_k != Ktot) {
        stop("Something went wrong: matrix of clusters has not been filled
             properly.")
    }

    # Save total number of elements of moc that are NA
    number_nas <- sum(is.na(moc))

    # If required
    if (fill) {
        if (fullData) {
            dataFill <- data
        } else {
            dataFill <- NULL
        }
        # Fill matrix of cluster labels by imputing missing labels based on the
        # available ones
        fillMOCoutput <-
            fillMOC(clLabels,
                    computeAccuracy = computeAccuracy,
                    data = dataFill)

        # Replace matrix of cluster labels
        clLabels <- fillMOCoutput$clLabels

        # Expand matrix of cluster labels and overwrite matrix-of-clusters
        moc <- expandMOC(fillMOCoutput$clLabels)
    }

    if (fill) {
        output <-
            list(
                moc = moc,
                datasetIndicator = datasetIndicator,
                number_nas = number_nas,
                clLabels = clLabels,
                K = K,
                fillMOCoutput = fillMOCoutput
            )
    } else {
        output <-
            list(
                moc = moc,
                datasetIndicator = datasetIndicator,
                number_nas = number_nas,
                clLabels = clLabels,
                K = K
            )
    }
    return(output)
}
