#' Plot silhouette
#'
#' Plot average silhouette for different values of K.
#' @param sil vector of the average silhouette for K from 2 to some value maxK
#' @param chosenK Chosen number of clusters. If specified, a vertical line is
#' plotted in correspondence of the indicated value. Default is NULL.
#' @param fileName Name of the png file.
#' @author Alessandra Cabassi \email{alessandra.cabassi@mrc-bsu.cam.ac.uk}
#' @keywords internal
#'
plotSilhouette <- function(sil, chosenK = NULL, fileName) {
    maxK <- length(sil) + 1

    fileName <- paste(fileName, "_silhouette.png", sep = "")
    grDevices::png(fileName, width = 400, height = 400)
    graphics::plot(2:maxK,
                   sil,
                   xlab = "Number of clusters",
                   ylab = "Average silhouette",
                   type = "o")
    if (!is.null(chosenK))
        graphics::abline(v = chosenK, col = "darkred")
    grDevices::dev.off()

    invisible(0)
}

#' Plot widest gap
#'
#' Plot widest gap for different values of K.
#' @param widestGap Vector of widest gap values for K from 2 to some value maxK.
#' @param chosenK Chosen number of clusters. If specified, a vertical line is
#' plotted in correspondence of the indicated value. Default is NULL.
#' @param fileName Name of the png file.
#' @author Alessandra Cabassi \email{alessandra.cabassi@mrc-bsu.cam.ac.uk}
#' @keywords internal
#'
plotWidestGap <- function(widestGap, chosenK = NULL, fileName) {
    maxK <- length(widestGap) + 1

    fileName <- paste(fileName, "_widestGap.png", sep = "")
    grDevices::png(fileName, width = 400, height = 400)
    graphics::plot(2:maxK,
                   widestGap,
                   xlab = "Number of clusters",
                   ylab = "Widest gap",
                   type = "o")
    if (!is.null(chosenK))
        graphics::abline(v = chosenK, col = "darkred")
    grDevices::dev.off()
}

#' Plot Dunn's index
#'
#' Plot Dunn's index (minimum separation / maximum diameter) for different
#' values of K.
#' @param dunns Vector of widest gap values for K from 2 to some value maxK.
#' @param chosenK Chosen number of clusters. If specified, a vertical line is
#' plotted in correspondence of the indicated value. Default is NULL.
#' @param fileName Name of the png file.
#' @author Alessandra Cabassi \email{alessandra.cabassi@mrc-bsu.cam.ac.uk}
#' @keywords internal
#'
plotDunns <- function(dunns, chosenK = NULL, fileName) {
    maxK <- length(dunns) + 1

    fileName <- paste(fileName, "_dunns.png", sep = "")
    grDevices::png(fileName, width = 400, height = 400)
    graphics::plot(2:maxK,
                   dunns,
                   xlab = "Number of clusters",
                   ylab = "Dunn's index",
                   type = "o")
    if (!is.null(chosenK))
        graphics::abline(v = chosenK, col = "darkred")
    grDevices::dev.off()

}

#' Plot Dunn's alternative index
#'
#' Plot Dunn's alternative index (minimum average dissimilarity between two
#' cluster divided by maximum average within cluster dissimilarity) for
#' different values of K.
#' @param dunns Vector of widest gap values for K from 2 to some value maxK.
#' @param chosenK Chosen number of clusters. If specified, a vertical line is
#' plotted in correspondence of the indicated value. Default is NULL.
#' @param fileName Name of the png file.
#' @author Alessandra Cabassi \email{alessandra.cabassi@mrc-bsu.cam.ac.uk}
#' @keywords internal
#'
plotDunn2s <- function(dunns, chosenK = NULL, fileName) {
    maxK <- length(dunns) + 1

    fileName <- paste(fileName, "_dunn2s.png", sep = "")
    grDevices::png(fileName, width = 400, height = 400)
    graphics::plot(2:maxK,
                   dunns,
                   xlab = "Number of clusters",
                   ylab = "Dunn's alternative index",
                   type = "o")
    if (!is.null(chosenK))
        graphics::abline(v = chosenK, col = "darkred")
    grDevices::dev.off()

    invisible(0)
}

#' Choose K that maximises the silhouette from a set of kernel matrices and
#' clusterings
#'
#' @param kernelMatrix N X N X (maxK-1) array of kernel matrices.
#' @param clLabels (maxK-1) X N matrix containing the clusterings obtained for
#' different values of K
#' @param maxK Maximum number of clusters considered.
#' @param savePNG If TRUE, a plot of the silhouette is saved in the working
#' folder. Defaults to FALSE.
#' @param fileName If savePNG is TRUE, this is the name of the png file.
#' @param isDistance Boolean. If TRUE, the kernel matrices are interpreted as
#' matrices of distances, otherwise as matrices of similarities.
#' @param widestGap Boolean. If TRUE, also computes widest gap index (and plots
#' it if savePNG is TRUE).
#' @param dunns Boolean. If TRUE, also computes Dunn's index: minimum separation
#' / maximum diameter (and plots it if savePNG is TRUE).
#' @param dunn2s Boolean. If TRUE, also computes an alternative version of
#' Dunn's index: minimum average dissimilarity between two cluster / maximum
#' average within cluster dissimilarity (and plots it if savePNG is TRUE).
#' @return The function returns a list containing `silh`, a vector of length
#' maxK-1 such that silh[i] is the silhouette for K = i+1, and `K`, the lowest
#' number of clusters for which the silhouette is maximised.
#' @author Alessandra Cabassi \email{alessandra.cabassi@mrc-bsu.cam.ac.uk}
#' @export
#'
maximiseSilhouette <-
    function(kernelMatrix,
             clLabels,
             maxK,
             savePNG = FALSE,
             fileName = "silhouette",
             isDistance = FALSE,
             widestGap = FALSE,
             dunns = FALSE,
             dunn2s = FALSE) {

    # Initialise vector of average silhouette
    sil <- rep(NA, maxK - 1)
    if (widestGap)
        widestGap_ <- rep(NA, maxK - 1)
    if (dunns)
        dunns_ <- rep(NA, maxK - 1)
    if (dunn2s)
        dunn2s_ <- rep(NA, maxK - 1)

    for (i in 2:maxK) {
        if (isDistance) {
            DM <-
                stats::as.dist(kernelMatrix[, , i - 1], diag = FALSE,
                               upper = FALSE)
        } else {
            # Use the kernel matrix as distance matrix
            DM <-
                stats::as.dist(1 - kernelMatrix[, , i - 1], diag = FALSE,
                               upper = FALSE)
        }

        # Calculate the average silhouette over all the points
        if (widestGap | dunns | dunn2s)
            cluster_stats <- fpc::cluster.stats(DM, clLabels[i - 1, ])
        if (widestGap)
            widestGap_[i - 1] <- cluster_stats$widestgap
        if (dunns)
            dunns_[i - 1] <- cluster_stats$dunn
        if (dunn2s)
            dunn2s_[i - 1] <- cluster_stats$dunn2
        sil[i - 1] <-
            mean(cluster::silhouette(clLabels[i - 1, ], DM)[, "sil_width"])

    }

    K <- which.max(sil)[1] + 1
    # The "+1" is there because the silhouette is stored in a vector where
    # element i corresponds to number of clusters equal to i+1
    if (savePNG) {
        plotSilhouette(sil, K, fileName)
        if (widestGap)
            plotWidestGap(widestGap_, K, fileName)
        if (dunns)
            plotDunns(dunns_, K, fileName)
        if (dunn2s)
            plotDunn2s(dunn2s_, K, fileName)

    }

    return(list(silhouette = sil, K = K))
}
