\name{marginTree.predict}
\alias{marginTree.predict}
\title{Prediction from the
 margin tree  classifier}
\description{
A function to predict from the
 margin tree  classifier
}
\usage{
marginTree.predict(train.obj,x, threshold=1)
}

\arguments{
\item{train.obj}{Output from call to marginTree}
  \item{x}{Data of feature values, n samples by p features}
  \item{threshold}{Threshold values for feature selection. Between 1 and 0.)}
}
\details{
  \code{marginTree.predict} yields class predictions from a margin tree classifier.}

\value{
Vector of predicted values
}

\references{
 Rob Tibshirani and Trevor Hastie. Tech report. Feb. 2006.
Margin trees for high-dimensional classification
  Available at http://www-stat.stanford.edu/~tibs/research.html} 


\author{Robert Tibshirani and Trevor Hastie}

\examples{
#generate some data with 5 classes and 100 features

set.seed(543)
x=matrix(rnorm(40*1000),nrow=40)
y=sort(rep(1:5,8))
x[y==2 | y==3, 1:50]=x[y==2|y==3, 1:50]+1
x[y==3,51:100]=x[y==3,51:100]+1
x[y==4|y==5,1:50]=x[y==4|y==5,1:50]-1
x[y==5, 51:100]=x[y==5,51:100]+1

#generate some test data
xx=x+5*matrix(rnorm(40*1000),nrow=40)
yy=y


#train the classifier

train.obj<-   marginTree(x,y)

# get predictions for test data
marginTree.predict(train.obj,xx)
}

