
#' Fang and Wang's instability measure
#'
#' This method is based on the algorithm developped by Fang and Wang but with more choice regarding the instability measure.
#' Their measure is equivalent to 1 - rand.index, here one can chose any normalized similarity measure and the instability will be 1 - similarity.
#' For each number of clusters, several pair of bootstrap subsambles are selected and the instability measure is computed from the clustering of these pairs.
#' The optimal number of clusters is the value for which the instability is the lowest.
#'
#' @param X data matrix or data frame of size n x d, n observations and d features
#' @param maxK maximum number of clusters to evaluate
#' @param B number of resampling iterations
#' @param clusterAlg clustering algorithm. Its output must be a list containing attributs "cluster" and "predict".
#' For more details, check the formatting of function \code{\link{myKmean}}.
#' @param similarity function measuring the similarity between two partitions.
#' @param verbose logical. If TRUE, plots the evolution of the algorithm
#' @param ... additional parameters for the clustering algorithm
#'
#' @return List with 3 components:
#' \describe{
#' \item{\preformatted{inst_mean}}{vector containing the mean instability measure for 2 to maxK cluster number}
#' \item{\preformatted{kopt}}{the optimal number of clusters}
#' \item{\preformatted{instability}}{matrix containing the instability measures for all cluster number and all subsampling iterations.}
#' }
#' @export
#' @importFrom fossil adj.rand.index
#'
#' @references  Fang, Y. and Wang, J. (2012). Selection of the number of clusters via the bootstrap method. Computational Statistics  Data Analysis, 56:468-477.
bootstrapInstability <- function(X, maxK, B = 50, clusterAlg = myKmean, similarity = adj.rand.index, verbose =TRUE,...){
  n <- dim(X)[1]
  instability <- matrix(0, maxK-1, B)
  cluster1 <- cluster2 <- numeric(n)
  for (b in 1:B){
    if(verbose==TRUE){
      print(paste("Iteration: ", b))
    }
    idx1 <- sample(1:n,n,replace = TRUE)
    idx2 <- sample(1:n,n,replace = TRUE)
    inv_idx1 <- setdiff(1:n,idx1)
    inv_idx2 <- setdiff(1:n,idx2)
    X1 <- X[idx1,]
    X2 <- X[idx2,]
    for (k in 2:maxK){
      X1.cluster <- clusterAlg(X1, k,...)
      X2.cluster <- clusterAlg(X2, k,...)
      cluster1[idx1] <- X1.cluster$cluster
      cluster2[idx2] <- X2.cluster$cluster
      cluster1[inv_idx1] <- X1.cluster$predict(X[inv_idx1,])
      cluster2[inv_idx2] <- X2.cluster$predict(X[inv_idx2,])
      instability[k-1,b] <- 1 - similarity(cluster1, cluster2)
    }
  }
  inst_mean <- apply(instability,1,mean)
  plot(2:maxK, inst_mean, main = "Instability", xlab = "number of clusters", ylab = "instability", type = "b")
  return(list("inst_mean"= inst_mean, "kopt"=which.min(inst_mean)+1, "instability"=instability))
}
