---
title: "Introduction to coca"
author: "Alessandra Cabassi"
date: "2019-06-03"
output:
  html_document:
    toc: true
    theme: united
vignette: >
  %\VignetteIndexEntry{Introduction to coca}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

# Introduction

The R package **coca** contains the functions needed to use COCA (**C**luster-**O**f-**C**lusters **A**nalysis), an integrative clustering method that was first introduced in a breast cancer study by The Cancer Genome Atlas in 2012 and quickly became a popular tool in cancer studies (see e.g. Hoadley et al. 2014 and Aure et al. 2017). It is based on \emph{Consensus Clustering} (Monti et al., 2013), an algorithm that was initially developed to assess the stability of clusters obtained with any clustering algorithm.

The main goal of COCA is to summarise clusterings found in different ’omics datasets, by identifying a “global" clustering across the datasets that is intended to be in good agreement with the clustering structures
identified in each of the individual datasets.

# Building the matrix of clusters

The first step of COCA is the construction of the Matrix-Of-Clusters (MOC). This is a binary matrix of size _N x K_, where _K_ is the sum of the number of clusters _Km_ in every dataset _Xm_. Therefore, to each column _j_ of this matrix corresponds a cluster _mk_ in dataset _Xm_. The _(i,j)_-th entry of the matrix of clusters is equal to one if data point _i_ belongs to cluster _mk_ in dataset _Xm_, and is equal to zero otherwise.For more details about how the matrix of clusters is built, please see Cabassi and Kirk (2018).

The function that can be used to build a matrix of clusters starting from a list of heterogeneous datasets (referring to the same observations) is `buildMOC`.
In the example below, we assume to have three datasets with the same number of clusters, six. The clustering structure in each of those dataasets is found via k-means clustering.


```r
### Load data
data <- list()
data[[1]] <- as.matrix(read.csv(system.file("extdata",
                      "dataset1.csv", package = "coca"), row.names = 1))
data[[2]] <- as.matrix(read.csv(system.file("extdata",
                      "dataset2.csv", package = "coca"), row.names = 1))
data[[3]] <- as.matrix(read.csv(system.file("extdata",
                      "dataset3.csv", package = "coca"), row.names = 1))

### Build matrix of clusters
outputBuildMOC <- coca::buildMOC(data, M = 3, K = 6, distances = "cor")

### Extract matrix of clusters and dataset indicator vector
moc <- outputBuildMOC$moc
datasetIndicator <- outputBuildMOC$datasetIndicator
```

The package also contains a function that can be used to plot the resulting Matrix-Of-Clusters, `plotMOC`. Here we use as annotations the true cluster labels, but in real applications annotations can be any dataframe with
one element for each row of the matrix of clusters (the row names must correspond!). Please be aware that each
column of the dataframe that contains categorical variables must be defined with `as.factor()` if you want
each category to have a different colour (otherwise they will be treaded as continuous variables and each
category will have a different shade of the same colour).


```r
### Prepare annotations
true_labels <- as.matrix(read.csv(system.file("extdata", "cluster_labels.csv",
                package = "coca"), row.names = 1))
annotations <- data.frame(true_labels = as.factor(true_labels))

### Plot matrix of clusters
coca::plotMOC(moc, datasetIndicator, annotations = annotations)
```

![plot of chunk plot_moc](figure/plot_moc-1.png)

Here the datasets don't have names, so they have been assigned integer numbers. If available, you can specify
cluster names; this will make the row names easier to interpret.


```r
### Prepare annotations
true_labels <- as.matrix(read.csv(system.file("extdata", "cluster_labels.csv",
                package = "coca"), row.names = 1))
annotations <- data.frame(true_labels = as.factor(true_labels))

### Set dataset names
datasetNames <- c(rep("A", 6), rep("B", 6), rep("C", 6))

### Plot matrix of clusters
coca::plotMOC(moc, datasetIndicator, datasetNames = datasetNames,
        annotations = annotations)
```

![plot of chunk plot_moc_with_names](figure/plot_moc_with_names-1.png)

As you can see, the first part of each row name corresponds to the dataset, the second one to the cluster index. Moreover, each colour in the main matrix corresponds to one dataset (here dataset A is green, dataset B is yellow, and so on).   

# Cluster-Of-Clusters Analysis

The MOC matrix is then used as input to consensus clustering (CC), an algorithm that was developed by Monti et al. (2003) to assess cluster stability when analysing a single dataset. The resulting consensus matrix is then used as the similarity matrix for a hierarchical clustering method (or any other distance-based clustering algorithm). These last two steps of the COCA algorithm are contained in the `coca()` function.


```r
### COCA

# Use COCA to find global clustering
coca <- coca(moc, K = 6, hclustMethod = "average")

# Compare clustering to the true labels
ari <- mclust::adjustedRandIndex(true_labels, coca$clusterLabels)
ari
```

```
## [1] 0.8163432
```

```r
### Plot the matrix of clusters with the newly found cluster labels
annotations$coca <- as.factor(coca$clusterLabels)
coca::plotMOC(moc, datasetIndicator, datasetNames = datasetNames,
        annotations = annotations)
```

![plot of chunk coca](figure/coca-1.png)

Again, if the number of clusters is not know a priori, the user can delegate the choice of K to the `coca()` function. There are two methods available here: the silhouette, where the distance between data points _i_ and _j_ is defined as _1_ minus the _(i,j)_-th kernel entry in the final kernel matrix, and the delta area under the curve, which is the method suggested by Monti et al. (2003). Please note that the properties of the latter
have not been assessed yet.


```r
# Use COCA to find global clustering and chooose the number of clusters
coca <- coca(moc, maxK = 10, hclustMethod = "average")

# Compare clustering to the true labels
ari <- mclust::adjustedRandIndex(true_labels, coca$clusterLabels)
ari
```

```
## [1] 0.8163432
```

```r
### Plot the matrix of clusters with the newly found cluster labels
annotations$coca <- as.factor(coca$clusterLabels)
coca::plotMOC(moc, datasetIndicator, datasetNames = datasetNames,
        annotations = annotations)
```

![plot of chunk coca_unknownK](figure/coca_unknownK-1.png)

# References

Aure, M.R., Vitelli, V., Jernström, S., Kumar, S., Krohn, M., Due, E.U., Haukaas, T.H., Leivonen, S.K., Vollan, H.K.M., Lüders, T. and Rødland, E. (2017). Integrative clustering reveals a novel split in the luminal A subtype of breast cancer with impact on outcome. Breast Cancer Research, 19(1), p.44.

Cabassi, A. and Kirk, P. D. W. (2019). Multiple kernel learning for integrative consensus clustering of genomic datasets. arXiv preprint. arXiv:1904.07701.

Hoadley, K.A., Yau, C., Wolf, D.M., Cherniack, A.D., Tamborero, D., Ng, S., Leiserson, M.D., Niu, B., McLellan, M.D., Uzunangelov, V. and Zhang, J., 2014. Multiplatform analysis of 12 cancer types reveals molecular classification within and across tissues of origin. Cell, 158(4), pp.929-944.

Monti, S. et al. (2003). Consensus Clustering: A Resampling-Based Method for Class Discovery and Visualization of Gene. Machine Learning, 52(i), 91–118.

The Cancer Genome Atlas (2012). Comprehensive molecular portraits of human breast tumours. Nature,
487(7407), 61–70.
