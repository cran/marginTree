\name{marginTree.plotcv}
\alias{marginTree.plotcv}
\title{A function to plot the cross-validated  error curves from
  classifier 
a margin tree classifier}
\description{A function to plot the cross-validated error curves
from a margin tree classifier}
\usage{
marginTree.plotcv(cv.obj)
}

\arguments{
  \item{cv.obj}{The result of a call to marginTree.cv}
}    

\details{
  \code{marginTree.plotcv} plots the cross-validated misclassification error
  curves,   from a margin tree classifier.
   An overall plot, and a plot by class, are produced.
}
\value{
}

\references{}


\author{Robert Tibshirani and Trevor Hastie}

\examples{

set.seed(543)

# generate some data 
x=matrix(rnorm(40*1000),nrow=40)
y=sort(rep(1:5,8))
x[y==2 | y==3, 1:50]=x[y==2|y==3, 1:50]+1
x[y==3,51:100]=x[y==3,51:100]+1
x[y==4|y==5,1:50]=x[y==4|y==5,1:50]-1
x[y==5, 51:100]=x[y==5,51:100]+1


#train the classifier

train.obj<-  marginTree(x,y)


# do cross-validation
cv.obj<- marginTree.cv(x,y,train.obj)


#make cv plot
marginTree.plotcv(cv.obj)
}

\keyword{ }

