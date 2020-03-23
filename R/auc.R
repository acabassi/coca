#' Compute Area Under the Curve (AUC)
#'
#' This function allows to compute the area under the curve of the
#' empirical distribution function of a consensus matrix
#' as described in Monti et al. (2003), Section 3.3.1.
#'
#' @param consensusMatrix Consensus matrix, output of the
#' "coca::consensusCluster" function.
#' @return This function returns the area under the curve.
#' @author Alessandra Cabassi \email{alessandra.cabassi@mrc-bsu.cam.ac.uk}
#' @keywords internal
#'

computeAUC <- function(consensusMatrix) {
    N <- dim(consensusMatrix)[1]
    x <- sort(as.vector(consensusMatrix))
    Z <- length(x)
    empirical_cdf <- stats::ecdf(x)
    auc <-
        sum((x[seq_len(Z - 1) + 1] - x[seq_len(Z - 1)]) *
                empirical_cdf(x[seq_len(Z - 1) + 1]))

}

#' Plot area under the curve
#'
#' Plot area under the curve for different values of K
#' @param deltaAUC Vector of the difference between the area under the curve
#' between each value K of the number of clusters and K-1. For the smallest
#' value considered (usually two) this is simply the area under the curve for
#' that value of cluster number.
#' @param chosenK Chosen number of clusters. If specified, a vertical line is
#' plotted in correspondance of the indicated value. Default is NULL.
#' @param fileName name of the png file
#' @return invisible(0)
#' @author Alessandra Cabassi \email{alessandra.cabassi@mrc-bsu.cam.ac.uk}
#' @keywords internal
#'
plotDeltaAUC <-
    function(deltaAUC,
             chosenK = NULL,
             fileName = "deltaAUC.png") {

    maxK <- length(deltaAUC) + 1

    if (!dir.exists("delta-auc"))
        dir.create("delta-auc", showWarnings = FALSE)

    fileName <- paste(fileName, ".png", sep = "")
    grDevices::png(fileName, width = 400, height = 400)
    graphics::plot(2:maxK,
                   deltaAUC,
                   xlab = "Number of clusters",
                   ylab = "Relative change in area under the curve",
                   type = "o")
    if (!is.null(chosenK))
        graphics::abline(v = chosenK)
    grDevices::dev.off()

    invisible(0)
}

#' Choose number of clusters based on AUC
#'
#' This function allows to choose the number of clusters in a dataset
#' based on the area under the curve of the empirical distribution
#' function of a consensus matrix, calculated for different (consecutive)
#' cluster numbers, as explained in the article by Monti et al. (2003),
#' Section 3.3.1.
#'
#' @param areaUnderTheCurve Vector of length maxK-1 containing the area
#' under the curve of the empirical distribution function of the
#' consensus matrices obtained with K varying from 2 to maxK.
#' @param savePNG Boolean. If TRUE, a plot of the area under the curve
#' for each value of K is saved as a png file. The file is saved in a
#' subdirectory of the working directory, called "delta-auc". Default is FALSE.
#' @param fileName If "savePNG" is TRUE, this is the name of the png file.
#' Can be used to specify the folder path too. Default is "deltaAUC". The ".png"
#' extension is automatically added to this string.
#' @return This function returns a list containing "deltaAUC", a vector of
#' length maxK-1 where element i is  the area under the curve for
#' K = i+1 minus the area under the curve for K = i (for i = 2 this
#' is simply the area under the curve for K = i), and "K" the lowest
#' among the values of K that are chosen by the algorithm.
#' @author Alessandra Cabassi \email{alessandra.cabassi@mrc-bsu.cam.ac.uk}
#' @examples
#' # Assuming that we want to choose among any value of K (number of clusters)
#' # between 2 and 10 and that the area under the curve is as follows:
#' areaUnderTheCurve <- c(0.05, 0.15, 0.4, 0.5, 0.55, 0.56, 0.57, 0.58, 0.59)
#'
#' # The optimal value of K can be chosen with:
#' K <- chooseKusingAUC(areaUnderTheCurve)$K
#' @references Monti, S., Tamayo, P., Mesirov, J. and Golub, T., 2003. Consensus
#' clustering: a resampling-based method for class discovery and visualization
#' of gene expression microarray data. Machine learning, 52(1-2), pp.91-118.
#' @export
#'
chooseKusingAUC <-
    function(areaUnderTheCurve,
             savePNG = FALSE,
             fileName =
                 "deltaAUC.png") {

    # Get value of maximum number of clusters considered
    maxK <- length(areaUnderTheCurve) + 1

    # Initialise vector of AUC[i]-AUC[i-1]
    deltaAUC <- rep(NA, maxK - 1)

    # For K=2, this cannot be computed so it is simply AUC for K = 2
    deltaAUC[1] <- areaUnderTheCurve[1]

    # Since the values in vector `areaUnderTheCurve` are not always
    # monotonically increasing, we need to store at each step the maximum value
    # encountered so far
    maxAUC <- areaUnderTheCurve[1]

    # Fill in vector deltaAUC according to Equation 7 in Monti et al. (2003)
    for (i in 2:(maxK - 1)) {
        deltaAUC[i] <- (areaUnderTheCurve[i] - maxAUC)/maxAUC
        maxAUC <- max(areaUnderTheCurve[i], maxAUC)
    }

    # Choose the value K such that deltaAUC[K+1] - deltaAUC[K] is smallest (not
    # its absolute value)
    K <- max(which(deltaAUC > 0.025)) + 1

    if (savePNG)
        plotDeltaAUC(deltaAUC, K, fileName)

    output <- list(deltaAUC = deltaAUC, K = K[1])
    return(output)
}
