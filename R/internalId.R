internalValid <- function(X, maxK, clusterAlg = myKmean, indexes = c("silhouette", "ch", "dunn") ){
  id <- 1:10
  names(id) <- c("silhouette", "ch", "dunn")
  indexAlg <- c(sil, ch, )
}


sil <- function(X, cluster){
  k <- length(table(cluster))
  score <- 0
  for (i in 1:k){
    Xi <- X[cluster==i,]
    ni <- nrow(Xi)
    a <- apply(as.matrix(dist(Xi)),2,sum) / (ni-1)
    b <- array(Inf, ni)
    for (j in (1:k)[-i]){
      Xj <- X[cluster==j,]
      nj <- nrow(Xj)
      tmp <- rowMeans(proxy::dist(Xi,Xj))
      b <- apply(rbind(b,tmp), 2, min)
    }
    score <- score + sum((b-a) / apply(rbind(a,b), 2, max))
  }
  score <- score / nrow(X)
  return(score)
}

ch <- function(X, cluster){
  n <- nrow(X)
  k <- length(table(cluster))
  nk <- table(cluster)
  d <- mean(dist(X)^2)
  dk <- aggregate(1:n, list(cluster), function(i) if(is.matrix(X[i,])){return(mean(dist(X[i,])^2))} else{return(0)})[-1]
  W <- 0.5*sum((nk-1)*dk)
  A <- sum((nk-1)*(d-dk))
  B <- 0.5*((k-1)*d +A)
  return(B * (n-k)/(W * (k-1)))
}

dunn <- function(X, cluster){

}
