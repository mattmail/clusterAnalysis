internalValid <- function(X, maxK, clusterAlg = myKmean, indexes = c("silhouette", "ch", "dunn") ){
  id <- 1:10
  names(id) <- c("silhouette", "ch", "dunn")
  indexAlg <- c(sil, ch, )
}


sil <- function(X, cluster){
  return(mean(silhouette(cluster, dist(X))[,3]))
}

ch <- function(X, cluster){
  n <- nrow(X)
  d <- ncol(X)
  nk <- table(cluster)
  wk <- aggregate(1:n, list(cluster), function(i) if(is.matrix(X[i,])){return(cov(X[i,]))} else{return(0)})[-1]
  W <- matrix(apply((nk-1)*wk, 2, sum), nrow = d)
  B <- (n-1)*cov(X) - W
  return(sum(diag(B)) * (n-k)/(sum(diag(W)) * (k-1)))
}

db <- function(X, cluster, p=2){
  n <- nrow(X)
  d <- ncol(X)
  k <- length(table(cluster))
  c <- aggregate(1:n, list(cluster), function(i) if(is.matrix(X[i,])){return(apply(X[i,], 2, mean))} else{return(X[i,])})[-1]
  dk <- matrix(0,k)
  c <- matrix(0, nrow = k, ncol = d)
  for(i in 1:k){
    idx <- cluster == i
    c[i,] <-  apply(X[idx, ], 2, mean)
    dk[i] <- (sum(abs(t(X[idx, ]) - as.numeric(c[i,]))^p)/ sum(idx))^(1/p)
  }
  m <- as.matrix(dist(c, p=p))
  diag(m) <- array(Inf, k)
  tf <- matrix(1, k)
  didj <- tf %*% t(dk) + dk %*% t(tf)
  db <- mean(apply(didj/m, 2, max))
  return(db)
}
