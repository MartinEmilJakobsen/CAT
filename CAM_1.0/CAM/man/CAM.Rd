\name{CAM}
\alias{CAM}
\title{
Causal Additive Model
}
\description{
fits a causal additive model using the CAM algorithm, see references below
}
\usage{
CAM(X, scoreName = "SEMGAM", parsScore = list(numBasisFcts = 10), numCores = 1, 
    maxNumParents = min(dim(X)[2] - 1, round(dim(X)[1]/20)), output = FALSE, 
    variableSel = FALSE, variableSelMethod = selGamBoost, 
    variableSelMethodPars = list(atLeastThatMuchSelected = 0.02, 
    atMostThatManyNeighbors = 10), pruning = FALSE, pruneMethod = selGam, 
    pruneMethodPars = list(cutOffPVal = 0.001, numBasisFcts = 10), intervData = FALSE, 
    intervMat = NA)
}
\arguments{
  \item{X}{
nxp matrix of training inputs (n data points, p dimensions)
}
  \item{scoreName}{
specifies the model type which is used to compute the score. Default is "SEMGAM" which assumes a generalized additive model class. Other options include "SEMLIN" which fits a linear model. 
}
  \item{parsScore}{
additional parameters can be supported to the score function.
}
  \item{numCores}{
specifies the number of cores that can be used for computation.
}
  \item{maxNumParents}{
specifies the maximal number of parents that are allowed in the model.
}
  \item{output}{
shall output be printed to the console (TRUE/FALSE)
}
  \item{variableSel}{
specifies whether initial variable selection (Step 1 of CAM algorithm) shall be performed (TRUE) or not (FALSE). Initial variable selection reduces the number of possible parents for a given node and therefore enables computing the causal structure for large p.
}
  \item{variableSelMethod}{
specifies the method that is used for variable selection. Default is selGamBoost which uses the gamboost function from mboost package. Other options include: selGam (gam() from mgcv), selLm based on linear regression, selLasso based on Lasso regression from package glmnet.
}
  \item{variableSelMethodPars}{
optional parameters to modify settings of the selection method.
}
  \item{pruning}{
specifies whether pruning (Step 3 of CAM algorithm) shall be performed (TRUE) or not (FALSE). Pruning reduces the number of edges in the estimated causal structure.
}
  \item{pruneMethod}{
specifies the method used for the pruning step. Default is selGAM which is based on the gam() function from the mgcv package.
}
  \item{pruneMethodPars}{
optional parameters to tune the pruning step.
}
  \item{intervData}{
boolean that indicates whether we use interventional data.
}
  \item{intervMat}{
the matrix intervMat has the same dimension as X. entry (i,j) == TRUE indicates that in experiment i, variable j has been intervened on. 
}
}
\details{
The code fits a CAM model. See the references below for more details. Identifiability results for the model class can be found in

J. Peters, J. Mooij, D. Janzing, B. Sch\"olkopf: Causal Discovery with Continuous Additive Noise Models, JMLR 15:2009-2053, 2014.
}
\value{
list of attributes of the final estimated causal structure
\item{Adj}{adjacency matrix of estimated causal graph}
\item{Score}{Total edge score of estimated graph}
\item{timesVec}{Vector containing various time measurements for execution times of the individual steps of the CAM algorithm}
}
\references{
P. B\"uhlmann, J. Peters, J. Ernest: 
CAM: Causal Additive Models, high-dimensional Order Search and Penalized Regression
Annals of Statistics 42:2526-2556, 2014.
}
\author{
Jonas Peters <jonas.peters@tuebingen.mpg.de> and Jan Ernest <ernest@stat.math.ethz.ch>
}


\examples{
n <- 500
eps1<-rnorm(n)
eps2<-rnorm(n)
eps3<-rnorm(n)
eps4<-rnorm(n)

x2 <- 0.5*eps2
x1 <- 0.9*sign(x2)*(abs(x2)^(0.5))+0.5*eps1
x3 <- 0.8*x2^2+0.5*eps3
x4 <- -0.9*sin(x3) - abs(x1) + 0.5*eps4

X <- cbind(x1,x2,x3,x4)

trueDAG <- cbind(c(0,1,0,0),c(0,0,0,0),c(0,1,0,0),c(1,0,1,0))
## x4 <- x3 <- x2 -> x1 -> x4
## adjacency matrix:
## 0 0 0 1
## 1 0 1 0
## 0 0 0 1
## 0 0 0 0

estDAG <- CAM(X, scoreName = "SEMGAM", numCores = 1, output = TRUE, variableSel = FALSE, 
              pruning = TRUE, pruneMethod = selGam, pruneMethodPars = list(cutOffPVal = 0.001))

cat("true DAG:\n")
show(trueDAG)

cat("estimated DAG:\n")
show(estDAG$Adj)
}

\keyword{Causality }
\keyword{Regression }
\keyword{Additive Noise Models }
\keyword{Restricted Structural Equation Model}

