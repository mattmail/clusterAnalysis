Gap <- function(X, maxK, clusterAlg = myKmean, B = 50, ...){
  X <- as.matrix(X)
  W <- numeric(maxK+1)
  n <- nrow(X)
  W[1] <- sum(dist(X))/n
  for (k in 2:(maxK+1)){
    cluster <- clusterAlg(X, k)$cluster
    idx <- 1:n
    W[k] <- sum(aggregate(idx, list(cluster), function(i)sum(dist(X[i,]))/length(i))[-1])
  }

  Wb <- matrix(0, maxK+1, B)
  for (b in 1:B){
    Xtmp <- scale(X, center = TRUE, scale = FALSE)
    v <- svd(Xtmp)$v
    Xtmp <- Xtmp %*% v
    rng <- apply(Xtmp, 2, range)
    Xb <- apply(rng, 2, function(r) runif(n=n, r[1], r[2]))
    Xb <- Xb %*% t(v)
    Wb[1,b] <- sum(dist(Xb))/n
    for (k in 2:(maxK+1)){
      cluster <- clusterAlg(Xb, k)$cluster
      idx <- 1:n
      Wb[k,b] <- sum(aggregate(idx, list(cluster), function(i)sum(dist(Xb[i,]))/length(i))[-1])
    }
  }
  lbar <- apply(log(Wb), 1, mean)
  gap <- lbar - log(W)
  sdk <- sqrt(apply((log(Wb)-lbar)^2, 1, mean))
  s <- sdk * sqrt(1+1/B)
  k <- 1
  while (k < maxK+1 && gap[k] < gap[k+1] - s[k+1]) {
    k <- k+1
  }
  if(k == maxK+1){
    kopt <- 1
  }
  else{
    kopt <- k
  }
  return(list(kopt=kopt, gap = gap, s = s))
  }
