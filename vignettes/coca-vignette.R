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

