#' Plot Matrix-Of-Clusters
#'
#' This function creates a matrix of clusters, starting from a list of heterogeneous
#' datasets.
#'
#' @param moc Matrix-Of-Clusters of size N x sumK.
#' @param datasetIndicator Vector containing
#' @param datasetNames Vector containing the names of the datasets to which each column
#' of labels corresponds. If NULL, datasetNames will be the same as datasetIndicator.
#' Default is NULL.
#' @param annotations Dataframe containing annotations. Number of rows must be N.
#' If the annotations are integers, use `as.factor()` for a better visual result.
#' @param clr Cluster rows. Default is FALSE.
#' @param clc Cluster columns. Default is FALSE.
#' @param save Boolean. If TRUE, plot is saved as a png file.
#' @param fileName File name for the plot if save is TRUE. Default is "moc.png".
#' @param showObsNames Boolean. If TRUE, the plot will also include the column names
#' (i.e. name of each observation). Default is FALSE, since there are usually too
#' many columns.
#' @author Alessandra Cabassi \email{ac2051@cam.ac.uk}
#' @references The Cancer Genome Atlas, 2012. Comprehensive molecular portraits of
#' human breast tumours. Nature, 487(7407), pp.61–70.
#' @examples
#' ## Load data
#' data <- list()
#' data[[1]] <- as.matrix(read.csv(system.file("extdata", "dataset1.csv",
#' package = "coca"), row.names = 1))
#' data[[2]] <- as.matrix(read.csv(system.file("extdata", "dataset2.csv",
#' package = "coca"), row.names = 1))
#' data[[3]] <- as.matrix(read.csv(system.file("extdata", "dataset3.csv",
#' package = "coca"), row.names = 1))
#'
#' ## Build matrix of clusters
#' outputBuildMOC <- buildMOC(data, M = 3, K = 6)
#'
#' ## Extract matrix of clusters and dataset indicator vector
#' moc <- outputBuildMOC$moc
#' datasetIndicator <- outputBuildMOC$datasetIndicator
#'
#' ## Prepare annotations
#' true_labels <- as.matrix(read.csv(system.file("extdata", "cluster_labels.csv",
#' package = "coca"), row.names = 1))
#' annotations <- data.frame(true_labels = as.factor(true_labels))
#'
#' ## Plot matrix of clusters
#' plotMOC(moc, datasetIndicator, annotations = annotations)
#'
#' @export
#'
plotMOC = function(moc, datasetIndicator, datasetNames = NULL, annotations = NULL,
                   clr = FALSE, clc = FALSE, save = FALSE, fileName = "moc.png",
                   showObsNames = FALSE){

    moc = t(moc)

    # Get number of datasets
    M = length(table(datasetNames))

    # Get sum of the number of clusters in each dataset
    sumK = dim(moc)[1]

    # If dataset names are not provided
    if(is.null(datasetNames)){

        # Dataset names = dataset indicators
        datasetNames = as.character(datasetIndicator)
    }

    M = length(table(datasetNames))

    # Make dataset names unique by adding a different number at the end of
    # the dataset name for each cluster
    datasetNamesLong <- c()
    for(i in 1:sumK){

        nClustersLeft <- sum(datasetNames[i:sumK] == datasetNames[i])
        datasetNamesLong <- c(datasetNamesLong, paste(datasetNames[i],
                              nClustersLeft, sep = " "))
    }

    # We want every dataset to have a different colour :)
    for(i in 1:sumK){
        moc[i,] <- moc[i,]*datasetIndicator[i]
    }

    # Plot!
    rownames(moc) <- datasetNamesLong
    # if(!is.null(annotations) && is.null(rownames(annotations))){
    #     rownames(annotations) <- colnames(moc)
    # }

    if(save) grDevices::png(fileName, width = 1000, height = 600)

    pheatmap::pheatmap(moc,  legend = TRUE,
           legend_breaks = 0:M,
           legend_labels = c("", unique(datasetNames)),
           color =  c("white", (RColorBrewer::brewer.pal(n = max(3,M), name = "Set3"))),
           cluster_rows = clr, clustering_distance_rows = "binary",
           cluster_cols = clc, clustering_distance_cols = "binary",
           annotation_col = annotations, show_colnames = showObsNames,
           drop_levels = FALSE, na_col = "seashell2")

    if(save){
        grDevices::dev.off()
        warning('After saving a pheatmap plot to file, you sometimes have to repeat the
                `dev.off()` command in order to shut down the plotting device completely.')
    }
}
