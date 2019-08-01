
#' Last Leap and Last Major Leap
#'
#' This method first computes the minimum distance between two cluster centers d.
#' Then, it gets the optimal number of clusters according to the last leap method and the last major leap method.
#'
#' @param X data matrix or data frame of size n x d, n observations and d features
#' @param maxK maximum number of cluster to evaluate
#' @param clusterAlg clustering algorithm. Its output must be a list having an attribute "centers" containing the centers of each cluster.
#'  For more details, check the formatting of function \code{\link{myKmean}}.
#' @param verbose logical, if TRUE it plot the evolution of the algorithm
#' @param ... additional parameters for the clustering algorithm
#'
#' @return list with 5 compoments:
#' \describe{
#' \item{\preformatted{d}}{vector of the minimum inter-center distance for each number of cluster}
#' \item{\preformatted{ll}}{vector containing the relative difference between d(k) and d(k+1)}
#' \item{\preformatted{ll_kopt}}{optimal number of clusters according to the last leap method}
#' \item{\preformatted{lml}}{vector indicating wether there is a major gap between values d(k) and d(k+1)}
#' \item{\preformatted{lml_kopt}}{optimal number of clusters according to the last major leap method}
#' }
#' @importFrom stats dist
#' @export
#' @references Gupta, A., Datta, S., and Das, S. (2018). Fast automatic estimation of the number of clusters from the minimum inter-center distancefor k-means clustering.Pattern Recognition Letters, 116.
#'
gupta <- function(X, maxK, clusterAlg = myKmean, verbose = FALSE, ...){
  d <- numeric(maxK)
  # cluster the data for each k and compute the minimum inter-center distance
  for (k in 2:(maxK+1)){
    if(verbose){
      print(paste("clustering with k =",k))
    }
    centers <- clusterAlg(X,k,...)$centers
    d[k-1] <- min(dist(centers))^2
  }
  ll <- (d[1:(maxK-1)]-d[2:maxK])/d[1:(maxK-1)]
  # shift ll so that the score for k clusters matches with kth position
  ll_tmp <- numeric(maxK)
  ll_tmp[2:maxK] <- ll
  ll <- ll_tmp
  ll_kopt <- which.max(ll)
  # Compute for each k maxd[k] the max of d[l] for l > k
  maxd <- d
  max <- d[maxK-1]
  for (i in (maxK-1):1){
    if(max < d[i+1]){
      max <- d[i+1]
    }
    maxd[i] <- max
  }
  lml <- 1*(0.5*d - maxd > 0)
  # shift lml so that the score for k clusters matches with kth position
  tmp <- 1:maxK
  tmp[2:maxK] <- lml[1:maxK-1]
  lml <- tmp
  # multiply lml with 1:maxK so that it is the last maximum that is returned
  lml_kopt <- which.max(1:maxK*lml)
  if(verbose){
    print(paste("According to ll the optimal number of clusters is", ll_kopt))
    print(paste("According to lml the optimal number of clusters is", lml_kopt))
  }
  return(list("d"=d,"ll"=ll, "lml"=lml, "ll_kopt" = ll_kopt, "lml_kopt" = lml_kopt))
}
