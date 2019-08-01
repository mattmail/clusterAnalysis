
#' Figure of Merit
#'
#' The figure of merit evaluate the quality of a clustering algorithm by computing it on all features but one and then, it calculates the mean squared error on the remaining feature for each cluster.
#' To find the optimal number of clusters, one has to find the "knee" in the curve.
#'
#' @param X data matrix or data frame of size n x d, n observations and d features
#' @param maxK maximum number of cluster to evaluate.
#' @param clusterAlg clustering algorithm. Its output must be a list having a compoment "cluster" containing the assignation of each observation.
#' For more details, check the formatting of function \code{\link{myKmean}}.
#' @param adjusted logical. If TRUE the adjusted figure of merit is returned.
#' @param verbose logical. If true, the figure of merit is plotted against the nuber of clusters.
#' @param ... additional parameters for the clustering algorithm.
#'
#' @return fom_scores, the list of scores
#' @importFrom stats aggregate
#' @export
#'
#' @references Yeung,  K.  Y.,  Haynor,  D.  R.,  and  Ruzzo,  W.  L.  (2001). Validating clustering for gene expression data. Bioinformatics,  17(4):309-318. \url{https://doi.org/10.1093/bioinformatics/17.4.309}
#'
fom <- function(X, maxK, clusterAlg = myKmean, adjusted = TRUE, verbose = TRUE, ...){
  n <- nrow(X)
  m <- ncol(X)
  fom_scores <- numeric(maxK)
  for (k in 1:maxK){
    fom_e_scores <- numeric(m)
    for (e in 1:m){
      clusters <- clusterAlg(X[ ,-e], k,...)$cluster
      fom_e_scores[e] <- sqrt(sum(aggregate(X[ ,e], list(clusters), function(x){sum((x - mean(x))^2)})[,-1])/n)
    }
    if (adjusted){
      fom_scores[k] <- sqrt(n/(n-k)) * sum(fom_e_scores)
    }
    else{
      fom_scores[k] <- sum(fom_e_scores)
    }
  }
  if(verbose){
    plot(fom_scores)
  }
  return(fom_scores)
}
