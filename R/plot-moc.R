#' Plot Matrix-Of-Clusters
#'
#' This function creates a matrix of clusters, starting from a list of
#' heterogeneous datasets.
#'
#' @param moc Matrix-Of-Clusters of size N x sumK.
#' @param datasetIndicator Vector containing integers indicating which rows
#' correspond to some clustering of the same dataset.
#' @param datasetNames Vector containing the names of the datasets to which each
#' column of labels corresponds. If NULL, datasetNames will be the same as
#' datasetIndicator. Default is NULL.
#' @param annotations Dataframe containing annotations. Number of rows must be
#' N. If the annotations are integers, use \code{as.factor()} for a better
#' visual result.
#' @param clr Cluster rows. Default is FALSE.
#' @param clc Cluster columns. Default is FALSE.
#' @param savePNG Boolean. If TRUE, plot is saved as a png file.
#' @param fileName If \code{savePNG} is TRUE, this is the string containing the
#' name of the moc figure. Can be used to specify the folder path too. Default
#' is "moc". The ".png" extension is automatically added to this string.
#' @param showObsNames Boolean. If TRUE, the plot will also include the column
#' names (i.e. name of each observation). Default is FALSE, since there are
#' usually too many columns.
#' @param showClusterNames Boolean. If TRUE, plot cluster names next to
#' corresponding row. Default is FALSE.
#' @param annotation_colors Optional. See annotation_colors in
#' pheatmap::pheatmap.
#' @author Alessandra Cabassi \email{alessandra.cabassi@mrc-bsu.cam.ac.uk}
#' @references The Cancer Genome Atlas, 2012. Comprehensive molecular portraits
#' of human breast tumours. Nature, 487(7407), pp.61–70.
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
#' # Create vector of dataset names, in the same order as they appear above
#' datasetNames <- c("Dataset1", "Dataset2", "Dataset3")
#'
#' # Build matrix of clusters
#' outputBuildMOC <- buildMOC(data, M = 3, K = 6, distances = "cor")
#'
#' # Extract matrix of clusters and dataset indicator vector
#' moc <- outputBuildMOC$moc
#' datasetIndicator <- outputBuildMOC$datasetIndicator
#'
#' # Prepare annotations
#' true_labels <- as.matrix(read.csv(
#' system.file("extdata", "cluster_labels.csv", package = "coca"),
#' row.names = 1))
#' annotations <- data.frame(true_labels = as.factor(true_labels))
#'
#' # Plot matrix of clusters
#' plotMOC(moc,
#'         datasetIndicator,
#'         datasetNames = datasetNames,
#'         annotations = annotations)
#' @export
#'
plotMOC <-
    function(moc,
             datasetIndicator,
             datasetNames = NULL,
             annotations = NULL,
             clr = FALSE,
             clc = FALSE,
             savePNG = FALSE,
             fileName = "moc.png",
             showObsNames = FALSE,
             showClusterNames = FALSE,
             annotation_colors = NA) {

    moc <- t(moc)

    # Get number of datasets
    M <- length(table(datasetNames))

    # Get sum of the number of clusters in each dataset
    sumK <- dim(moc)[1]

    # If dataset names are not provided
    if (is.null(datasetNames)) {
        # Dataset names = dataset indicators
        datasetNames <- as.character(datasetIndicator)
    }

    M <- length(table(datasetNames))

    # Make dataset names unique by adding a different number at the end of the
    # dataset name for each cluster
    datasetNamesLong <- c()
    for (i in seq_len(sumK)) {
        nClustersLeft <- sum(datasetNames[i:sumK] == datasetNames[i])
        datasetNamesLong <-
            c(datasetNamesLong, paste(datasetNames[i], nClustersLeft,
                                      sep = " "))
    }

    # We want every dataset to have a different colour :)
    for (i in seq_len(sumK)) {
        moc[i, ] <- moc[i, ] * datasetIndicator[i]
    }

    # Plot!
    rownames(moc) <- datasetNamesLong
    # if(!is.null(annotations) && is.null(rownames(annotations))){
    # rownames(annotations) <- colnames(moc) }


    if (M == 2) {
        mycols <- c("white", (RColorBrewer::brewer.pal(n = 2, name = "RdBu")))
        mycols <- mycols[c(1, 2, 4)]
    } else {
        mycols <- c("white", (RColorBrewer::brewer.pal(n = max(3, M),
                                                       name = "Set3")))
    }

    if (savePNG)
        grDevices::png(paste(fileName,".png"), width = 1000, height = 600)

    pheatmap::pheatmap(moc,
                       legend = TRUE,
                       legend_breaks = 0:M,
                       legend_labels = c("", unique(datasetNames)),
                       color = mycols,
                       cluster_rows = clr,
                       clustering_distance_rows = "binary",
                       cluster_cols = clc,
                       clustering_distance_cols = "binary",
                       annotation_col = annotations,
                       show_colnames = showObsNames,
                       annotation_colors = annotation_colors,
                       show_rownames = showClusterNames,
                       drop_levels = FALSE,
                       na_col = "seashell2",
                       border_color = NA)

    if (savePNG)
        grDevices::dev.off()
}
