#' Plot Matrix-Of-Clusters
#'
#' This function creates a matrix of clusters, starting from a list of heterogeneous datasets.
#'
#' @param moc Matrix-Of-Clusters of size N x sumK.
#' @param datasetIndicator Vector containing
#' @param datasetNames Vector containing the names of the datasets to which each column of
#' labels corresponds. If NULL, datasetNames will be the same as datasetIndicator.
#' Default is NULL.
#' @param annotations Vector or matrix of annotations of size N x a where a is the
#' number of different annotations.
#' @param clr Cluster rows. Default is FALSE.
#' @param clc Cluster columns. Default is FALSE.
#' @param save Boolean. If TRUE, plot is saved as a pdf file.
#' @param fileName File name. Default is "moc.pdf".
#' @author Alessandra Cabassi \email{ac2051@cam.ac.uk}
#' @references The Cancer Genome Atlas, 2012. Comprehensive molecular portraits of human breast
#' tumours. Nature, 487(7407), pp.61–70.
#' @export
#'
plotMOC = function(moc, datasetIndicator, datasetNames, annotations = NA,
                   clr = FALSE, clc = FALSE, save = FALSE, fileName = "moc.pdf"){

  moc = t(moc)

  # Get number of datasets
  M = length(table(datasetNames))

  # Get sum of the number of clusters in each dataset
  sumK = dim(moc)[1]

  # If dataset names are not provided
  if(is.null(datasetNames)){
    # Dataset names = dataset indicators
    datasetNames = as.string(datasetIndicator)
  }

  # Make dataset names unique by adding a different number at the end of
  # the dataset name for each cluster
  for(i in 1:sumK){
    nClustersLeft <- sum(datasetNames == datasetNames[i])
      datasetNames[i] = paste(datasetNames[i], nClustersLeft, sep = " ")
  }

  # We want every dataset to have a different colour :)
  for(i in 1:sumK){
    moc[i,] <- moc[i,]*datasetIndicator[i]
  }


  # Plot!
  rownames(moc) <- datasetNames
  if(save) pdf(fileName, width = 10, height = 5)
  pheatmap::pheatmap(moc,  legend = TRUE,
           legend_breaks = 0:M,
          # legend_labels = c("0", table(datasetNames)),
           color =  c("white", (RColorBrewer::brewer.pal(n = max(3,M), name = "Set3"))),
           cluster_rows = clr, clustering_distance_rows = "binary",
           cluster_cols = clc, clustering_distance_cols = "binary",
           annotation_col = annotations)

  if(save){
    dev.off()
    warning('After saving a pheatmap plot to file, you sometimes have to repeat the `dev.off()` command
            in order to shut down the plotting device completely.')
  }

}
