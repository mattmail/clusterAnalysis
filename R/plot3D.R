#' 3D dataset visualization
#'
#' Function for plotting the 3 principal components of a dataset. One can also add color labels to visualise cluster assignations.
#'
#' @param X data matrix or data frame of size n x d, n observations and d features. d must be superior or equal to 3.
#' @param label array of size n containing each cluster assignation for each observation.
#'
#' @export
#' @importFrom grDevices rainbow
#' @importFrom car scatter3d
#'
plot3D <- function(X, label = NULL){
  X <- as.matrix(X)
  if(ncol(X) < 3){
    stop("The number of features must be more than 2")
  }
  X.pca <- prcomp(X, scale=F)$x
  if(is.null(label)){
    scatter3d(X.pca[,1],X.pca[,2],X.pca[,3],  point.col = "blue", surface=FALSE)
  }
  else{
    k <- length(table(label))
    scatter3d(X.pca[,1],X.pca[,2],X.pca[,3], point.col =rainbow(k)[label],surface=FALSE, grid = FALSE)
  }
}
