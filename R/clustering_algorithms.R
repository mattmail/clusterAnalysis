#' K-means clustering
#'
#' This function is a wrapper to \code{\link[stats:kmeans]{kmeans}} in the package stats.
#'
#' @param X data matrix or data frame of size n x d, n observations and d features
#' @param k number of clusters
#'
#' @importFrom stats kmeans
#' @return list of 3 components:
#' \describe{
#' \item{\preformatted{cluster}}{vector of integer between 1 and k containing the allocation of each point}
#' \item{\preformatted{centers}}{matrix (d x k) of the centers of each cluster}
#' \item{\preformatted{predict}}{function predicting to which cluster an observation belongs. The input can be a single observation vector or a matrix of several observations. The assigned cluster is the one with the closest center to the given observation.}
#' }
#' @export
#'
myKmean <- function(X, k){
  res <- kmeans(X,k, iter.max = 50)
  list("cluster" = res$cluster, "predict" = function(Z){
    closest.point <- function(x) {
    cluster.dist <- apply(res$centers, 1, function(y) sqrt(sum((x-y)^2)))
    return(which.min(cluster.dist)[1])}
    return(apply(Z, 1, closest.point))
  }, "centroids" = res$centers)
}


#' Spectral clustering
#'
#' This function is a wrapper to \code{\link[anocva:spectralClustering]{spectralClustering}} in the package anocva.
#'
#' @param X data matrix or data frame of size n x d, n observations and d features
#' @param k number of clusters
#' @param simi similarity measure between two vectors
#'
#' @importFrom anocva spectralClustering
#' @importFrom class knn
#' @return list of 3 components:
#' \describe{
#' \item{\preformatted{cluster}}{vector of integer between 1 and k containing the allocation of each point}
#' \item{\preformatted{centers}}{matrix (d x k) of the centers of each cluster}
#' \item{\preformatted{predict}}{function predicting to which cluster an observation belongs. The prediciton is done with the k-nearest neighbours function.}
#' }
#' @export
#'
mySpectralClustering <- function(X, k, simi){
  A <- simi(X)
  cluster <- spectralClustering(A,k)
  d <- dim(X)[2]
  centers <- aggregate(X, rep(list(cluster),d), mean)[,-(1:d)]
  prediction <- function(x){
    knn(X,x,cluster,k=10)
  }
  return(list("cluster" = cluster, "predict" = prediction, "centers"=centers))
}


#' Hierarchical clustering
#'
#' This function performs hierarchical clustering in either average, complete or single linkage.
#'
#' @param X data matrix or data frame of size n x d, n observations and d features
#' @param k number of clusters
#' @param distance character for the type of distance. Default is "euclidean"
#' @param linkage "single", "average" or "complete"
#'
#' @import stats
#' @importFrom class knn
#' @return list of 3 components:
#' \describe{
#' \item{\preformatted{cluster}}{vector of integer between 1 and k containing the allocation of each point}
#' \item{\preformatted{centers}}{matrix (d x k) of the centers of each cluster}
#' \item{\preformatted{predict}}{function predicting to which cluster an observation belongs. The prediciton is done with the k-nearest neighbours function.}
#' }
#' @export
#'
myHierarchicalClustering <- function(X,k, distance = "euclidean", linkage = "average"){
  dist_M <- dist(X, method = distance)
  h_clust <- hclust(dist_M, method = linkage)
  cluster <- cutree(h_clust, k=k)
  d <- dim(X)[2]
  centers <- aggregate(X, rep(list(cluster),d), mean)[,-(1:d)]

  prediction <- function(x){
    knn(X,x,cluster,k=10)
  }
  return(list("cluster" = cluster, "predict" = prediction, "centers"=centers))
}

#' Partitioning around medoids
#'
#' This function is a wrapper to \code{\link[cluster:pam]{pam}} in the package cluster.
#'
#' @param X data matrix or data frame of size n x d, n observations and d features
#' @param k number of clusters
#' @param ... additional parameter for pam function
#'
#' @importFrom cluster pam
#' @importFrom class knn
#' @return list of 3 components:
#' \describe{
#' \item{\preformatted{cluster}}{vector of integer between 1 and k containing the allocation of each point}
#' \item{\preformatted{centers}}{matrix (d x k) of the centers of each cluster}
#' \item{\preformatted{predict}}{function predicting to which cluster an observation belongs. The prediciton is done with the k-nearest neighbours function.}
#' }
#' @export
#'
myPam <- function(X,k, ...){
  pa <- pam(X,k,...)
  prediction <- function(x){
    knn(X,x,pa$clustering, k=10)
  }
  return(list("cluster"=pa$clustering, "centers"=pa$medoids, "predict"= prediction))
}
