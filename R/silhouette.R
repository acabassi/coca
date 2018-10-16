#' Plot silhouette
#'
#' Plot average silhouette.
#' @param sil vector of the average silhouette for K from 2 to some value maxK
#' @param chosenK Chosen number of clusters. If specified, a vertical line is plotted in
#' correspondance of the indicated value. Default is NULL.
#' @param fileName name of the png file
#' @author Alessandra Cabassi \email{ac2051@cam.ac.uk}
plotSilhouette = function(sil, chosenK = NULL, fileName){

    maxK = length(sil) + 1

    if(!dir.exists("silhouette")) dir.create("silhouette", showWarnings = FALSE)

    fileName = paste('silhouette/', fileName, '.png', sep = '')
    grDevices::png(fileName, width = 400, height = 400)
    graphics::plot(2:maxK, sil, xlab = "Number of clusters", ylab = "Average silhouette", type = "o")
    if(!is.null(chosenK)) graphics::abline(v = chosenK)
    grDevices::dev.off()

}

#' Choose K that maximises the silhouette from a set of kernel matrices and clusterings
#'
#' @param kernelMatrix N X N X (maxK-1) array of kernel matrices.
#' @param clLabels (maxK-1) X N matrix containing the clusterings obtained for different values of K
#' @param maxK Maximum number of clusters considered.
#' @param savePNG If TRUE, a plot of the silhouette is saved in the working folder. Defaults to FALSE.
#' @param fileName If savePNG is TRUE, this is the name of the png file.
#' @return The function returns a list containing `silh`, a vector of length maxK-1 such that
#' silh[i] is the silhouette for K = i+1, and `K`, the lowest number of clusters for which the
#' silhouette is maximised.
#' @author Alessandra Cabassi \email{ac2051@cam.ac.uk}
#' @examples
#' data <- as.matrix(read.csv(system.file("extdata", "dataset1.csv", package = "klic"), row.names = 1))
#'
#' cm_2cl <- consensusCluster(data, 2)
#' cm_3cl <- consensusCluster(data, 3)
#'
#' km <- array(NA, c(300, 300, 2))
#' km[,,1] <- klic::spectrumShift(cm_2cl, coeff = 1.1)
#' km[,,2] <- klic::spectrumShift(cm_3cl, coeff = 1.1)
#'
#' clLabels <- array(NA, c(2,300))
#'
#' # Use kernel k-means to divide data into two clusters
#' parameters_kkmeans <- list()
#' parameters_kkmeans$cluster_count <- 2
#' kkm <- klic::kkmeans(km[,,1], parameters_kkmeans) # km[,,1] because 2 is the first number we try
#' clLabels[1,] <- kkm$clustering
#'
#' # Use kernel k-means to divide data into three clusters
#' parameters_kkmeans$cluster_count <- 3
#' kkm <- klic::kkmeans(km[,,2], parameters_kkmeans) # km[,,2] because 3 is the second number we try
#'clLabels[2,] <- kkm$clustering
#'
#' # Call maximiseSilhouette function
#' maxSil <- maximiseSilhouette(km, clLabels, maxK = 3)
#' # The output of maximiseSilhouette contains:
#' # * the values of the average silhoeutte for each number of clusters
#' print(maxSil$silhouette)
#' # * the number of clusters that maximises the average silhouette
#' # (if more than one value for the number of clusters k maximises the silhouette,
#' # all the such values are reported in the output)
#' print(maxSil$k)
#' @export
maximiseSilhouette = function(kernelMatrix, clLabels, maxK,
                              savePNG = FALSE, fileName = "silhouette"){

  # Initialise vector of average silhouette
  sil <- rep(5, maxK-1)

  for(i in 2:maxK){

    # Use the kernel matrix as distance matrix
    DM <- stats::as.dist(1 - kernelMatrix[,,i-1], diag = FALSE, upper = FALSE)

    # Calculate the average silhouette over all the points
    sil[i-1]<- mean(cluster::silhouette(clLabels[i-1,], DM)[,'sil_width'])
  }

  K = which.max(sil)[1]+1
  # The "+1" is there because the silhouette is stored in a vector where element i corresponds to number of clusters equal to i+1
  if(savePNG) plotSilhouette(sil, K, fileName)

  return(list(silhouette = sil,  K = K))
}
