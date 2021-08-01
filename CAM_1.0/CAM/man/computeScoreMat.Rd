\name{computeScoreMat}
\alias{computeScoreMat}
\title{
auxiliary file for CAM: Computes the initial score matrix.
}
\description{
auxiliary file to CAM. Computes the initial score matrix.
}
\usage{
computeScoreMat(X, scoreName, numParents, output, numCores, selMat, parsScore, 
                intervMat, intervData)
}
\arguments{
  \item{X}{
nxp matrix of training inputs (n data points, p dimensions)
}
  \item{scoreName}{
specifies the model type which is used to compute the score. Default is "SEMGAM" which assumes a generalized additive model class. Other options include "SEMLIN" which fits a linear model. 
}
  \item{numParents}{
indicates how many parents we consider. If numParents = 1 (default), then the score matrix is of dimension (p-1) x p. If numParents = 2, then the score matrix is of dimension (p-1)(p-2) x p and so on
}
  \item{output}{
boolean indicating whether information about the progress is written to the console.
}
  \item{numCores}{
specifies the number of cores that can be used for computation.
}
  \item{selMat}{
indicating the possible parent relationships.
}
  \item{parsScore}{
additional parameters can be supported to the score function.
}
  \item{intervMat}{
the matrix intervMat has the same dimension as X. entry (i,j) == TRUE indicates that in experiment i, variable j has been intervened on. 
}
  \item{intervData}{
boolean that indicates whether we use interventional data.
}
}
\value{
A list with elements
\item{scoreMat}{The score matrix. scoreMat[i,j] contains the gain in score if we consider i being a parent of j }
\item{rowParents}{Contains the row names of the score matrix. Only relevant if numParents > 1.}
\item{scoreEmtpyNodes}{Vector containing the scores of each node in the empty graph without any edges.}
}
\references{
P. B\"uhlmann, J. Peters, J. Ernest: 
CAM: Causal Additive Models, high-dimensional Order Search and Penalized Regression
Annals of Statistics 42:2526-2556, 2014.
}
\author{
J. Peters (jonas.peters@tuebingen.mpg.de) and J. Ernest (ernest@stat.math.ethz.ch)
}
\note{
This is an auxiliary file for CAM.
}

\seealso{
\code{\link[CAM]{CAM}}
}
