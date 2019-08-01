#' Clest: A prediction-based resampling method for estimating the number of clusters in a dataset
#'
#' This method measures the similarity of two clustering computed on two non-overlapping subsamples of the dataset.
#' The obtained measures are then compared with other similarity measures taken on generated datasets.
#' Those generated datasets are computed as in the gap statistics method (Tibshirani et al, 2001).
#'
#'
#' @param X data matrix or data frame of size n x d, n observations and d features
#' @param maxK maximum number of clusters to evaluate.
#' @param clusterAlg clustering algorithm. Its output must be a list having a compoment "cluster" containing the assignation of each observation.
#' For more details, check the formatting of function \code{\link{myKmean}}.
#' @param similarity function measuring the similarity between two partitions.
#' @param pmax threshold for the p-value
#' @param dmin threshold for the d-value
#' @param B number of resampling iterations
#' @param B0 number of reference datasets to generate
#' @param rho proportion of the train/test set
#' @param verbose logical, if TRUE, plots the evolution of the algorithm
#' @param ... additional parameters for the clustering algorithm
#'
#' @return list of 3 attributes:
#' \describe{
#' \item{\preformatted{p}}{vector of p-values}
#' \item{\preformatted{d}}{vector of d-values}
#' \item{\preformatted{kopt}}{optimal number of clusters}
#' }
#' @importFrom fossil adj.rand.index
#' @export
#'
#' @references
#' \itemize{
#' \item{Dudoit, S. and Fridlyand, J. (2002). A prediction-based resampling method for estimating the number of clusters in a dataset. Genome Biology, 3(7):research0036.1.
#' \url{https://doi.org/10.1186/gb-2002-3-7-research0036}}
#' \item{Tibshirani, R., Walther, G., and Hastie, T. (2001). Estimating the number of clusters in a data set via the gap statistic.Journal of the Royal Statistical Society Series B, 63:411-423.}}
clest <- function(X, maxK, clusterAlg = myKmean, similarity = adj.rand.index, pmax =0.05, dmin = 0.05, B = 50, B0 = 20, rho = 0.6, verbose = TRUE, ...){
  n <- nrow(X)
  X <- as.matrix(X)
  Skb <- matrix(0, maxK-1,B)
  if(verbose){
    print("Original dataset")
  }
  for (k in 2:maxK){
    if (verbose){
      print(paste("Number of clusters: ", k))
    }
    for (b in 1:B){
      idx_train <- sort(sample(1:n, floor(rho*n), replace = F))
      idx_test <- setdiff(1:n, idx_train)
      Xtrain <- X[idx_train,]
      Xtest <- X[idx_test,]
      cltrain <- clusterAlg(Xtrain,k,...)
      cltest <- clusterAlg(Xtest,k,...)
      testPred <- cltrain$predict(Xtest)
      Skb[k-1,b] <- similarity(cltest$cluster, testPred)
    }
  }
  t <- apply(Skb, 1, median)
  #repeat the same steps for the B0 reference datasets
  Xtmp <- scale(X, center = TRUE, scale = FALSE)
  v <- svd(Xtmp)$v
  Xtmp <- Xtmp %*% v
  rng <- apply(Xtmp, 2, range)
  tkb <- matrix(0, maxK-1, B0)
  for (b0 in 1:B0){
    if(verbose){
      cat('\n')
      print(paste("Reference dataset number", b0))
    }
    Zref <- apply(rng, 2, function(r, nn) runif(nn, min = r[1], max = r[2]), nn = n)
    Xref <- Zref %*% t(v)
    Skb_ref <- matrix(0, maxK-1,B)
    for (k in 2:maxK){
      if (verbose){
        print(paste("Number of clusters: ", k))
      }
      for (b in 1:B){
        idx_train <- sort(sample(1:n, floor(rho*n), replace = F))
        idx_test <- setdiff(1:n, idx_train)
        Xtrain <- Xref[idx_train,]
        Xtest <- Xref[idx_test,]
        cltrain <- clusterAlg(Xtrain,k,...)
        cltest <- clusterAlg(Xtest,k,...)
        testPred <- cltrain$predict(Xtest)
        Skb_ref[k-1,b] <- similarity(cltest$cluster, testPred)
      }
    }
    tkb[,b0] <- apply(Skb_ref, 1, median)
  }
  t0 <- apply(tkb, 1, mean)
  p <- apply(tkb >= t, 1, sum)/B0
  d <- t - t0
  optimalSet <- d[p <= pmax & d >= dmin]
  if(identical(optimalSet, numeric(0))){
    kopt <- 1
  }
  else{
    maxd <- max(optimalSet)
    kopt <- match(maxd, d) +1
  }
  return(list("p"=p, "d"=d, "kopt"= kopt))
}
