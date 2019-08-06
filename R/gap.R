#' Gap statistics
#'
#' Tibshirani's gap statistic for the determination of the number of clusters. It computes the within cluster dispertion of the partition
#' and it compares it with the within cluster dispertion of generated datasets having similar statistics to the original.
#'
#' The within-cluster dispertion is the normalized sum for each cluster of the sum of the distance between each pair in a cluster.
#'
#' @param X data matrix or data frame of size n x d, n observations and d features
#' @param maxK maximum number of clusters to evaluate.
#' @param clusterAlg clustering algorithm. Its output must be a list having a compoment "cluster" containing the assignation of each observation.
#' For more details, check the formatting of function \code{\link{myKmean}}.
#' @param B number of reference datasets to generate
#' @param null_distrib type of the null hypothesis. Can either be "gaussian", "uniform" or "uniformity".
#' "gaussian" draws observations from a mulidimensional normal distribution with the same mean and variance as in the original dataset for each feature .
#' "uniform" draws uniformely observations in the range of each feature. "uniformity" draws observation from a uniform distribution as in gap statistics (Tibshirani et al. 2001).#' @param verbose logical, if TRUE it plots the evolution of the algorithm
#' @param ... additional parameters for the clustering algorithm
#'
#' @return
#' @export
#'
#' @examples
Gap <- function(X, maxK, clusterAlg = myKmean, B = 50, null_hypothesis = "gaussian", verbose = TRUE, ...){
  X <- as.matrix(X)
  W <- numeric(maxK+1)
  n <- nrow(X)
  W[1] <- sum(dist(X))/n
  if(verbose){
    print("Original dataset")
  }
  for (k in 2:(maxK+1)){
    cluster <- clusterAlg(X, k, ...)$cluster
    idx <- 1:n
    W[k] <- sum(aggregate(idx, list(cluster), function(i)sum(dist(X[i,]))/length(i))[-1])
  }

  Wb <- matrix(0, maxK+1, B)
  switch(null_hypothesis,
    gaussian = {
      param <- rbind(apply(X, 2, mean), apply(X,2,var))
      distrib <- mvrnorm
    },
    uniformity = {
      Xtmp <- scale(X, center = TRUE, scale = FALSE)
      v <- svd(Xtmp)$v
      Xtmp <- Xtmp %*% v
      param <- apply(Xtmp, 2, range)
      distrib <- runif
    },
    uniform = {
      param <- apply(X, 2, range)
      distrib <- runif
    })
  for (b in 1:B){
    if(verbose){
      print(paste("Reference dataset number: ", b))
    }

    Xb <- apply(param, 2, function(p) distrib(n=n, p[1], p[2]))
    if (null_hypothesis == "uniformity"){
      Xb <- Xb %*% t(v)
    }

    Wb[1,b] <- sum(dist(Xb))/n
    for (k in 2:(maxK+1)){
      cluster <- clusterAlg(Xb, k, ...)$cluster
      idx <- 1:n
      Wb[k,b] <- sum(aggregate(idx, list(cluster), function(i)sum(dist(Xb[i,]))/length(i))[-1])
    }
  }
  lbar <- apply(log(Wb), 1, mean)
  gap <- lbar - log(W)
  sdk <- sqrt(apply((log(Wb)-lbar)^2, 1, mean))
  s <- sdk * sqrt(1+1/B)
  kopt <- which(gap[-(maxK+1)] >= gap[-1]-s[-1])[1]
  if(length(kopt)==0){
    kopt <- 1
  }
  plot(gap, main = "Gap statistics with error bars", xlab = "number of clusters", ylab = "Gap", type = "b")
  arrows(1:(maxK+1), gap+s, 1:(maxK+1), gap-s, angle = 90, code = 3, length = 1/8)
  return(list(kopt=kopt, gap = gap, s = s))
  }
