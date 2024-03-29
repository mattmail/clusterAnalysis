% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fom.R
\name{fom}
\alias{fom}
\title{Figure of Merit}
\usage{
fom(X, maxK, clusterAlg = myKmean, adjusted = TRUE, verbose = TRUE,
  ...)
}
\arguments{
\item{X}{data matrix or data frame of size n x d, n observations and d features}

\item{maxK}{maximum number of cluster to evaluate.}

\item{clusterAlg}{clustering algorithm. Its output must be a list having a compoment "cluster" containing the assignation of each observation.
For more details, check the formatting of function \code{\link{myKmean}}.}

\item{adjusted}{logical. If TRUE the adjusted figure of merit is returned.}

\item{verbose}{logical. If true, the figure of merit is plotted against the nuber of clusters.}

\item{...}{additional parameters for the clustering algorithm.}
}
\value{
fom_scores, the list of scores
}
\description{
The figure of merit evaluate the quality of a clustering algorithm by computing it on all features but one and then, it calculates the mean squared error on the remaining feature for each cluster.
To find the optimal number of clusters, one has to find the "knee" in the curve.
}
\references{
Yeung,  K.  Y.,  Haynor,  D.  R.,  and  Ruzzo,  W.  L.  (2001). Validating clustering for gene expression data. Bioinformatics,  17(4):309-318. \url{https://doi.org/10.1093/bioinformatics/17.4.309}
}
