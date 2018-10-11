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
#' @param fileName File name. Default is "moc.pdf".
#' @author Alessandra Cabassi \email{ac2051@cam.ac.uk}
#' @references The Cancer Genome Atlas, 2012. Comprehensive molecular portraits of human breast
#' tumours. Nature,
#' 487(7407), pp.61–70.
#' @export
#'
plotMOC = function(moc, datasetIndicator, datasetNames, annotations = NA,
                   clr = FALSE, clc = FALSE, fileName = "moc.pdf"){

  moc = t(moc)
  M = length(table(datasetNames))
  sumK = dim(moc)[1]
  if(is.null(datasetNames)){
    newDatasetNames = as.string(datasetIndicator)
  }else{
    for(i in 1:sumK){
      nClustersLeft <- sum(datasetNames == datasetNames[i])
      if(nClustersLeft>1){
        datasetNames[i] = paste(datasetNames[i], nClustersLeft, sep = " ")
      }
    }
  }


  for(i in 1:sumK){
    moc[i,] <- moc[i,]*datasetIndicator[i]
  }

  rownames(moc) <- datasetNames

  pdf(fileName, width = 10, height = 5)
  pheatmap::pheatmap(moc,  legend = TRUE,
           legend_breaks = 0:M,
          # legend_labels = c("0", table(datasetNames)),
           color =  c("white", (RColorBrewer::brewer.pal(n = max(3,M), name = "Set3"))),
           cluster_rows = clr, clustering_distance_rows = "binary",
           cluster_cols = clc, clustering_distance_cols = "binary",
           annotation_col = annotations)

  dev.off()

}
