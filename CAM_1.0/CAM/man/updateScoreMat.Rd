\name{updateScoreMat}
\alias{updateScoreMat}
\title{
auxiliary file for CAM: Updates the score matrix after having added edge i -> j to the graph.
}
\description{
auxiliary file for CAM. Updates the score matrix after having added edge i -> j to the graph.
}
\usage{
updateScoreMat(scoreMat, X, scoreName, i, j, scoreNodes, Adj, output, numCores, 
               maxNumParents, parsScore, intervMat, intervData)
}
\arguments{
  \item{scoreMat}{
the current score matrix that has to be updated
}
  \item{X}{
nxp matrix of training inputs (n data points, p dimensions)
}
  \item{scoreName}{
specifies the model type which is used to compute the score. Default is "SEMGAM" which assumes a generalized additive model class. Other options include "SEMLIN" which fits a linear model. 
}
  \item{i}{
starting point of the edge i->j that has been added
}
  \item{j}{
end point of the edge i->j that has been added
}
  \item{scoreNodes}{
vector containing current scores of each node
}
  \item{Adj}{
adjacency matrix of the graph
}
  \item{output}{
boolean indicating whether information about the progress is written to the console.
}
  \item{numCores}{
specifies the number of cores that can be used for computation.
}
  \item{maxNumParents}{
specifies the maximal number of parents that are allowed in the model.
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
\item{scoreMat}{the updated score matrix.}
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
