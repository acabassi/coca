## ----preparation_coca, fig.show='hold', message=FALSE, warning=FALSE, cache=TRUE----
### Load data
data <- list()
data[[1]] <- as.matrix(read.csv(system.file("extdata", 
                      "dataset1.csv", package = "coca"), row.names = 1))
data[[2]] <- as.matrix(read.csv(system.file("extdata", 
                      "dataset2.csv", package = "coca"), row.names = 1))
data[[3]] <- as.matrix(read.csv(system.file("extdata", 
                      "dataset3.csv", package = "coca"), row.names = 1))

### Set constants
n_datasets <- 3
n_clusters <- 6
N <- 300
true_labels <- as.matrix(read.csv(system.file("extdata",
                         "cluster_labels.csv", package = "coca"), row.names = 1))

### Fill label matrix with clusterings found with the k-means clustering algorithm
labelMatrix <- array(NA, c(n_clusters, N, n_datasets))
for(i in 1:n_datasets){
  output <- kmeans(data[[i]], n_clusters)
  for(k in 1:n_clusters){
    labelMatrix[k,,i] <- (output$cluster==k)
  }
}

# Convert label matrix from logic to numeric matrix
labelMatrix <- labelMatrix*1

# Build MOC matrix
MOC <- rbind(labelMatrix[,,1],labelMatrix[,,2], labelMatrix[,,3])

## ----coca, fig.show='hold', message=FALSE, warning=FALSE, cache=TRUE-----
### COCA

# Use COCA to find global clustering
coca <- coca(t(MOC), K = 6, hclustMethod = "average")

# Compare clustering to the true labels
ari <- mclust::adjustedRandIndex(true_labels, coca$clusterLabels)
ari

