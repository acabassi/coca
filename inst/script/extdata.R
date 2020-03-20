n_elements_per_cluster <- 50
n_clusters <- 5

true_labels <- c(rep(1, n_elements_per_cluster),
                 rep(2, n_elements_per_cluster),
                 rep(3, n_elements_per_cluster),
                 rep(4, n_elements_per_cluster),
                 rep(5, n_elements_per_cluster))

dataset1 <- t(mvtnorm::rmvnorm(2, true_labels))
dataset2 <- t(mvtnorm::rmvnorm(2, true_labels))
dataset3 <- t(mvtnorm::rmvnorm(2, true_labels))

write.csv(dataset1, "inst/extdata/dataset1.csv")
write.csv(dataset2, "inst/extdata/dataset2.csv")
write.csv(dataset3, "inst/extdata/dataset3.csv")
write.csv(true_labels, "inst/extdata/cluster_labels.csv")
