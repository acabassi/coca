## ----build_coca, fig.show='hold', message=FALSE, warning=FALSE, cache=TRUE----
### Load data
data <- list()
data[[1]] <- as.matrix(read.csv(system.file("extdata",
                      "dataset1.csv", package = "coca"), row.names = 1))
data[[2]] <- as.matrix(read.csv(system.file("extdata",
                      "dataset2.csv", package = "coca"), row.names = 1))
data[[3]] <- as.matrix(read.csv(system.file("extdata",
                      "dataset3.csv", package = "coca"), row.names = 1))

### Build matrix of clusters
outputBuildMOC <- coca::buildMOC(data, M = 3, K = 5, distances = "cor")

### Extract matrix of clusters and dataset indicator vector
moc <- outputBuildMOC$moc
datasetIndicator <- outputBuildMOC$datasetIndicator

## ----plot_moc, fig.show='hold', message=FALSE, warning=FALSE, cache=TRUE------

### Prepare annotations
true_labels <- as.matrix(read.csv(system.file("extdata", "cluster_labels.csv",
                package = "coca"), row.names = 1))
annotations <- data.frame(true_labels = as.factor(true_labels))

### Plot matrix of clusters
coca::plotMOC(moc, datasetIndicator, annotations = annotations)

## ----plotmoc_wnames, fig.show='hold', message=FALSE, warning=FALSE, cache=TRUE----

### Prepare annotations
true_labels <- as.matrix(read.csv(system.file("extdata", "cluster_labels.csv",
                package = "coca"), row.names = 1))
annotations <- data.frame(true_labels = as.factor(true_labels))

### Set dataset names
datasetNames <- c(rep("A", 5), rep("B", 5), rep("C", 5))

### Plot matrix of clusters
coca::plotMOC(moc, datasetIndicator, datasetNames = datasetNames,
        annotations = annotations)

## ----coca, fig.show='hold', message=FALSE, warning=FALSE, cache=TRUE----------
### COCA

# Use COCA to find global clustering
coca <- coca::coca(moc, K = 5)

# Compare clustering to the true labels
ari <- mclust::adjustedRandIndex(true_labels, coca$clusterLabels)
ari

### Plot the matrix of clusters with the newly found cluster labels
annotations$coca <- as.factor(coca$clusterLabels)
coca::plotMOC(moc, datasetIndicator, datasetNames = datasetNames,
        annotations = annotations)

## ----coca_unknownK, fig.show='hold', message=FALSE, warning=FALSE, cache=TRUE----

# Use COCA to find global clustering and chooose the number of clusters
coca <- coca::coca(moc, maxK = 10, hclustMethod = "average")

# Compare clustering to the true labels
ari <- mclust::adjustedRandIndex(true_labels, coca$clusterLabels)
ari

### Plot the matrix of clusters with the newly found cluster labels
annotations$coca <- as.factor(coca$clusterLabels)
coca::plotMOC(moc, datasetIndicator, datasetNames = datasetNames,
        annotations = annotations)

