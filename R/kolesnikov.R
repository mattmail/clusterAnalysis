
#' Quantatization error modeling
#'
#'This method is based on the within cluster variance W. The idea is to estimate a parameter \emph{a} and to find the minimum of this function : \eqn{W(k)k^a}{W(k)k^a}.
#'
#' @param X data matrix or data frame of size n x d, n observations and d features
#' @param maxK maximum number of clusters to evaluate
#' @param clusterAlg clustering algorithm. Its output must be a list containing parameters "cluster" and "center".
#' For more details, check the formatting of function \code{\link{myKmean}}.
#'
#' @return list having 3 attributes:
#' \describe{
#' \item{\preformatted{kopt}}{optimal number of clusters}
#' \item{\preformatted{PCF}}{vector of score}
#' \item{\preformatted{a}}{value of the exponent \emph{a}}
#' }
#' @export
#'
#' @references Kolesnikov, A., Trichina, E., and Kauranne, T. (2015). Estimating the number of clusters in a numerical data set via quantizationerror modeling.Pattern Recognition, 48:941-952.
kolesnikov <- function(X, maxK, clusterAlg = myKmean, ...){
  X <- as.matrix(X)
  W <- array(0, maxK)
  for (k in 1:maxK){
    clustering <- clusterAlg(X, k, ...)
    centers <- clustering$centers
    dist <- X - centers[clustering$cluster,]
    W[k] <- sum(dist^2)
  }
  numerator <- maxK*sum(log(1:maxK)^2)-sum(log(1:maxK))^2
  denominator <- maxK*sum(log(1:maxK)*log(W)) - sum(log(1:maxK))*sum(log(W))
  a <- -2*numerator/denominator
  pcf <- W*(1:maxK)^(2/a)
  kopt <- which.min(pcf)
  return(list("kopt"= kopt,"PCF" = pcf, "a"=a))
}
