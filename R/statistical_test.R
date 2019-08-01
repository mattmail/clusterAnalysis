#' Statistical test for clustering relevance
#'
#' Statistical test for checking the clustering of a specified dataset is relevant. Several datasets are generated under a null hypothesis
#' and their distribution of nearest neighbours distances are compared with the one of the original dataset.
#'
#'
#' @param X data matrix or data frame of size n x d, n observations and d features
#' @param s number of reference datasets to generate
#' @param null_distrib type of the null hypothesis. Can either be "gaussian", "uniform" or "uniformity".
#' "gaussian" draws observations from a mulidimensional normal distribution with the same mean and variance as in the original dataset for each feature .
#' "uniform" draws uniformely observations in the range of each feature. "uniformity" draws observation from a uniform distribution as in gap statistics (Tibshirani et al. 2001).
#'
#' @return list of 2 components
#' \describe{
#' \item{\preformatted{U}}{vector containing the discrepancy measures.}
#' \item{\preformatted{pvalue}}{proportion of discrepancy measure of the generated datasets that are at least as large as the discrepancy measure of the original dataset.}
#' }
#' @export
#'
#' @importFrom spatstat nndist
#' @importFrom MASS mvrnorm
#' @importFrom stats runif
#' @importFrom stats prcomp
#' @importFrom stats var
#' @importFrom stats ecdf
#'
#' @references
#' \itemize{
#' \item{McShane,  L.  M.,  Radmacher,  M.  D.,  Freidlin,  B.,  Yu, R., Li, M.-C., and Simon, R. (2002). Methods for assessing reproducibility of clustering patterns observed in analyses of microarray data.Bioinformatics, 18(11):1462-1469. \url{https://doi.org/10.1093/bioinformatics/18.11.1462}}
#' \item{Tibshirani, R., Walther, G., and Hastie, T. (2001). Estimating the number of clusters in a data set via the gap statistic.Journal of the Royal Statistical Society Series B, 63:411-423.}}
#'
statistical_test <- function(X, s, null_distrib = "gaussian"){
  X <- as.matrix(X)
  Xpca <- prcomp(X)$x
  n <- nrow(Xpca)
  if(dim(X)[2]>2){
    Xpca <- Xpca[,1:3]
    d <- 3
  }
  else{
    d <- 2
  }
  nndistances <- nndist(Xpca, k=1)
  ecdf_nn <- ecdf(nndistances)
  min_dist <- min(nndistances)
  max_dist <- max(nndistances)

  if (null_distrib == "gaussian"){
    mu <- apply(Xpca, 2, mean)
    sig2 <- apply(Xpca, 2, var)
    param <- rbind(mu, sig2)
    distrib <- mvrnorm
  }
  else if(null_distrib == "uniform"){
    param <- apply(Xpca, 2, range)
    distrib <- runif
  }
  else if(null_distrib == "uniformity"){
    Xtmp <- scale(Xpca, center = TRUE, scale = FALSE)
    v <- svd(Xtmp)$v
    Xtmp <- Xtmp %*% v
    param <- apply(Xtmp, 2, range)
    distrib <- runif
  }
  else {
    stop("null_distrib must either be gaussian or uniform")
  }

  G <- list(ecdf_nn)
  for (i in 2:s){
    Z <- apply(param, 2, function(p) distrib(n=n, param[1], param[2]))
    if(null_distrib == "uniformity"){
      Z <- Z %*% t(v)
    }
    G[[i]] <- ecdf(nndist(Z, k=1))
  }

  estimate_G <- function(i, x){
    f_list <- G[-i]
    res <- 0
    for(f in f_list){
      res <- res + f(x)
    }
    res/(s-1)
  }
  U <- numeric(s)
  for (i in 1:s){
    x <- seq(min_dist, max_dist, length = 100)
    U[i] <- sum((G[[i]](x)-estimate_G(i, x))^2)*(x[2]-x[1])
  }
  plot(function(x){x}, col='blue')
  lines(estimate_G(1,sort(nndistances)), ecdf_nn(sort(nndistances)), col="red")
  return(list(U=U, pvalue = sum(U[-1]>=U[1])/(s-1) ))
}
