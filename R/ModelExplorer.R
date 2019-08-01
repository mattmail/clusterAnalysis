#' Stability measure for finding the number of clusters
#'
#' Model Explorer measures the stability of a given clustering method on subsambles of the original dataset.
#' It plots the cumulative distribution of the stability measure. The more the distribution is concentrated on the right, the better.
#'
#' @param X data matrix or data frame of size n x d, n observations and d features
#' @param maxK maximum number of clusters to evaluate.
#' @param similarity function measuring the similarity between two partitions.
#' @param clusterAlg clustering algorithm. Its output must be a list having a compoment "cluster" containing the assignation of each observation.
#' For more details, check the formatting of function \code{\link{myKmean}}.
#' @param rho numeric between 0 and 1. Proportional size of the subsamples.
#' @param B Number of resampling iterations.
#' @param verbose logical. If TRUE, it plot the evolution of the algorithm.
#' @param ... additional parameters for the clustering algorithm.
#'
#' @return Matrix of size maxK-1 x B containing the stability measure for each iteration.
#'
#' @export
#' @importFrom fossil adj.rand.index
#' @importFrom stats ecdf
#' @import graphics
#' @references Ben-Hur, A., Elisseeff, A., and Guyon, I. (2002). A stability based method for discovering structure in clustered data.Pacific Symposium on Biocomputing. Pacific Symposium on Biocomputing, 2002:6-17
#'
ModelExplorer <- function(X, maxK, similarity=adj.rand.index, clusterAlg = myKmean, rho = 0.8, B = 100, verbose = FALSE,  ...){
  n <- nrow(X)
  Skb <- matrix(0, maxK-1,B)
  for (b in 1:B){
    if (verbose){
      print(paste("Iteration: ", b))
    }
    for (k in 2:maxK){
      #sample in proportion rho observations of the dataset
      idx1 <- sort(sample(1:n, floor(rho*n), replace = F))
      idx2 <- sort(sample(1:n, floor(rho*n), replace = F))
      X1 <- X[idx1,]
      X2 <- X[idx2,]
      cl1 <- clusterAlg(X1,k,...)$cluster
      cl2 <- clusterAlg(X2,k,...)$cluster
      # compute the similarity on the intersection of the two subsets
      inter <- intersect(idx1,idx2)
      Skb[k-1,b] <- similarity(cl1[is.element(idx1, inter)], cl2[is.element(idx2, inter)])
    }
  }
  cols <- c("blue", "green", "orange", "magenta", "cyan", "red", "yellow", "gray","black","purple","pink","brown","blue4","coral","coral2","coral4","azure","darkblue")
  plot((1:1000/1000), ecdf(Skb[1,])((1:1000/1000)), type="l", col = cols[1], xlab = "similarity", ylab = "cumulative" )
  for (k in 2:(maxK-1)){
    lines((1:1000/1000), ecdf(Skb[k,])((1:1000/1000)), col= cols[k])
  }
  legend("topleft", as.character(2:maxK), fill=cols[1:10])
  return(Skb)
}



