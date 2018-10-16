#' Compute Area Under the Curve (AUC)
#'
#' This function allows to compute the area under the curve of the
#' empirical distribution function of a consensus matrix
#' as described in Monti et al. (2003), Section 3.3.1.
#'
#' @param consensusMatrix Consensus matrix, output of the `coca::consensusCluster`
#'        function
#' @return This function returns the area under the curve.
#' @author Alessandra Cabassi \email{ac2051@cam.ac.uk}
#' @references Monti, S., Tamayo, P., Mesirov, J. and Golub, T., 2003. Consensus clustering:
#' a resampling-based method for class discovery and visualization of gene expression microarray
#' data. Machine learning, 52(1-2), pp.91-118.
#'
computeAUC = function(consensusMatrix){

    N <- dim(consensusMatrix)[1]
    x <- sort(as.vector(consensusMatrix))
    Z <- length(x)
    empirical_cdf <- ecdf(x)
    auc <- sum((x[2:Z]-x[1:Z-1])*empirical_cdf(x[2:Z]))

}

#' Plot area under the curve
#'
#' Plot area under the curve for different values of K
#' @param areaUnderTheCurve vector of the average silhouette for K from 2 to some value maxK
#' @param chosenK Chosen number of clusters. If specified, a vertical line is plotted in
#' correspondance of the indicated value. Default is NULL.
#' @param fileName name of the png file
#' @author Alessandra Cabassi \email{ac2051@cam.ac.uk}
plotDeltaAUC = function(deltaAUC, chosenK = NULL, fileName){

    maxK = length(deltaAUC) + 1

    if(!dir.exists("delta-auc")) dir.create("delta-auc", showWarnings = FALSE)

    fileName = paste('delta-auc/', fileName, '.png', sep = '')
    grDevices::png(fileName, width = 400, height = 400)
    graphics::plot(2:maxK, deltaAUC, xlab = "Number of clusters", ylab = "Area under the curve", type = "o")
    if(!is.null(chosenK)) graphics::abline(v = chosenK)
    grDevices::dev.off()

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
#'        under the curve of the empirical distribution function of the
#'        consensus matrices obtained with K varying from 2 to maxK.
#' @return This function returns a list containing `deltaAUC`, a vector of
#'         length maxK-1 where element i is  the area under the curve for
#'         K = i+1 minus the area under the curve for K = i (for i = 2 this
#'         is simply the area under the curve for K = i), and `K` the lowest
#'         among the values of K that are chosen by the algorithm.
#' @export
#'
chooseKusingAUC = function(areaUnderTheCurve, savePNG, fileName){

    # Get value of maximum number of clusters considered
    maxK <- length(areaUnderTheCurve) + 1

    # Initialise vector of AUC[i]-AUC[i-1]
    deltaAUC <- rep(NA, maxK-1)

    # For K=2, this cannot be computed so it is simply AUC for K = 2
    deltaAUC[1] <- areaUnderTheCurve[1]

    # Since the values in vector `areaUnderTheCurve` are not always
    # monotonically increasing, we need to store at each step the
    # maximum value encountered so far
    maxAUC <- areaUnderTheCurve[1]

    #Fill in vector deltaAUC according to Equation 7 in Monti et al. (2003)
    for(i in 2:(maxK-1)){
        deltaAUC[i] <- areaUnderTheCurve[i]-maxAUC
        maxAUC <- max(areaUnderTheCurve[i], maxAUC)
    }

    # Choose the value K such that deltaAUC[K+1] - deltaAUC[K]
    # is smallest (not its absolute value)
    K <- which.min(deltaAUC[2:(maxK-1)]-deltaAUC[1:(maxK-2)]) +1

    if(savePNG) plotDeltaAUC(deltaAUC, K, fileName)

    output <- list(deltaAUC = deltaAUC, K = K[1])
    return(output)
}
