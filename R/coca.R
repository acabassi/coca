#' Cluster-Of-Clusters Analysis
#'
#' This function allows to do Cluster-Of-Clusters Analysis on a binary matrix
#' where each column is a clustering of the data, each row corresponds to a data
#' point and the element in position (i,j) is equal to 1 if data point i belongs
#' to cluster j, 0 otherwise.
#'
#' @param moc N X C data matrix, where C is the total number of clusters
#' considered.
#' @param K Number of clusters.
#' @param maxK Maximum number of clusters considered for the final clustering if
#' K is not known. Default is 6.
#' @param B Number of iterations of the Consensus Clustering step.
#' @param pItem Proportion of items sampled at each iteration of the Consensus
#' Cluster step.
#' @param hclustMethod Agglomeration method to be used by the hclust function to
#' perform hierarchical clustering on the consensus matrix. Can be "single",
#' "complete", "average", etc. For more details please see ?stats::hclust.
#' @param choiceKmethod Method used to choose the number of clusters if K is
#' NULL, can be either "AUC" (area under the curve, work in progress) or
#' "silhouette". Default is "silhouette".
#' @param ccClMethod Clustering method to be used by the Consensus Clustering
#' algorithm (CC). Can be either "kmeans" for k-means clustering or "hclust" for
#' hiearchical clustering. Default is "kmeans".
#' @param ccDistHC Distance to be used by the hiearchical clustering algorithm
#' inside CC. Can be "pearson" (for 1 - Pearson correlation), "spearman" (for 1-
#' Spearman correlation), or any of the distances provided in stats::dist()
#' (i.e. "euclidean", "maximum", "manhattan", "canberra", "binary" or
#' "minkowski"). Default is "euclidean".
#' @param maxIterKM Number of iterations for the k-means clustering algorithm.
#' Default is 1000.
#' @param verbose Boolean.
#' @param savePNG Boolean. Save plots as PNG files. Default is FALSE.
#' @param fileName If \code{savePNG} is TRUE, this is the string containing (the
#' first part of) the name of the output files. Can be used to specify the
#' folder path too. Default is "coca". The ".png" extension is automatically
#' added to this string.
#' @param widestGap Boolean. If TRUE, compute also widest gap index to choose
#' best number of clusters. Default is FALSE.
#' @param dunns Boolean. If TRUE, compute also Dunn's index to choose best
#' number of clusters. Default is FALSE.
#' @param dunn2s Boolean. If TRUE, compute also alternative Dunn's index to
#' choose best number of clusters. Default is FALSE.
#' @param returnAllMatrices Boolean. If TRUE, return consensus matrices for all
#' considered values of K. Default is FALSE.
#' @return This function returns a list containing:
#' \item{consensusMatrix}{a symmetric matrix where the element in position
#' (i,j) corresponds to the proportion of times that items i and j have been
#' clustered together and a vector of cluster labels.}
#' \item{clusterLabels}{the final cluster labels.}
#' \item{K}{the final number of clusters. If provided by the user, this is
#' the same as the input. Otherwise, this is the number of clusters selected via
#' the requested method (see argument \code{choiceKmethod}).}
#' \item{consensusMatrices}{if returnAllMatrices = TRUE, this array
#' also returned, containing the consensus matrices obtained for each of the
#' numbers of clusters considered by the algorithm.}
#' @author Alessandra Cabassi \email{alessandra.cabassi@mrc-bsu.cam.ac.uk}
#' @references The Cancer Genome Atlas, 2012. Comprehensive molecular
#' portraits of human breast tumours. Nature, 487(7407), pp.61–70.
#' @references Cabassi, A. and Kirk, P. D. W. (2019). Multiple kernel learning
#' for integrative consensus clustering of 'omic datasets. arXiv preprint.
#' arXiv:1904.07701.
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
#' outputBuildMOC <- buildMOC(data, M = 3, K = 5, distances = "cor")
#'
#' # Extract matrix of clusters
#' moc <- outputBuildMOC$moc
#'
#' # Do Cluster-Of-Clusters Analysis
#' outputCOCA <- coca(moc, K = 5)
#'
#' # Extract cluster labels
#' clusterLabels <- outputCOCA$clusterLabels
#'
#' @export

coca <-
    function(moc,
             K = NULL,
             maxK = 6,
             B = 1000,
             pItem = 0.8,
             hclustMethod = "average",
             choiceKmethod = "silhouette",
             ccClMethod = "kmeans",
             ccDistHC = "euclidean",
             maxIterKM = 1000,
             savePNG = FALSE,
             fileName = "coca",
             verbose = FALSE,
             widestGap = FALSE,
             dunns = FALSE,
             dunn2s = FALSE,
             returnAllMatrices = FALSE) {

    # Intialise output list
    output <- list()

    N <- dim(moc)[1]

    if (is.null(K) & choiceKmethod == "silhouette") {
        consensusMatrix <- array(NA, c(N, N, maxK - 1))
        clLabels <- array(NA, c(maxK - 1, N))

        for (i in seq_len(maxK - 1) + 1) {
            ### Step 1. Compute the consensus matrix
            consensusMatrix[, , i - 1] <-
                consensusCluster(moc, i, B, pItem, clMethod = ccClMethod,
                                 dist = ccDistHC)
            ### Step 2. Use hierarchical clustering on the consensus matrix
            distances <- stats::as.dist(1 - consensusMatrix[, , i - 1])
            hClustering <- stats::hclust(distances, method = hclustMethod)
            clLabels[i - 1, ] <- stats::cutree(hClustering, i)
        }

        K <- maximiseSilhouette(consensusMatrix, clLabels, maxK, savePNG,
                                fileName, widestGap = widestGap, dunns = dunns,
                                dunn2s = dunn2s)$K

    } else if (is.null(K) & choiceKmethod == "AUC") {
        consensusMatrix <- array(NA, c(N, N, maxK - 1))
        areaUnderTheCurve <- rep(NA, maxK - 1)

        for (i in seq_len(maxK - 1) + 1) {
            ### Step 1. Compute the consensus matrix ###
            consensusMatrix[, , i - 1] <-
                consensusCluster(moc, i, B, pItem, clMethod = ccClMethod,
                dist = ccDistHC, maxIterKM = maxIterKM)
            ### Step 2. Compute area under the curve ###
            areaUnderTheCurve[i - 1] <- computeAUC(consensusMatrix[, , i - 1])
        }

        K <- chooseKusingAUC(areaUnderTheCurve, savePNG, fileName)$K

    } else if (is.null(K)) {
        stop("Method to choose number of clusters has not been recognised.
            Please make sure that it is either `silhouette` or `AUC`.")
    } else {
        consensusMatrix <- NULL
    }

    if (verbose)
        print(paste("K =", K, sep = " "))

    ### Step 1. Compute the consensus matrix ###
    if (!is.null(consensusMatrix)) {
        output$consensusMatrix <- consensusMatrix[, , K - 1]
    } else {
        output$consensusMatrix <- consensusCluster(moc, K, B, pItem,
                                                   clMethod = ccClMethod,
                                                   dist = ccDistHC)
    }

    ### Step 2. Use hierarchical clustering on the consensus matrix ###
    distances <- stats::as.dist(1 - output$consensusMatrix)
    hClustering <- stats::hclust(distances, method = hclustMethod)
    output$clusterLabels <- stats::cutree(hClustering, K)

    output$K <- K

    if (returnAllMatrices)
        output$consensusMatrices <- consensusMatrix

    return(output)

}
