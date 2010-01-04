\name{marginTree.cv}
\alias{marginTree.cv}
\title{Cross-validation of the
 margin tree  classifier}
\description{
A function to  cross-validate the
 margin tree  classifier}
\usage{
marginTree.cv(x, y, train.obj, nfold = min(table(y)), folds = NULL, threshold = NULL, n.threshold=20)
}

\arguments{
  \item{x}{Data of feature values, n samples by p features}
  \item{y}{Vector of class labels}
  \item{train.obj}{Output from call to \code{marginTree}}
  \item{nfold}{Number of cross-validation folds (default 10)}
\item{folds}{A list with nfold components, each component a vector of
    indices of the samples in that fold. By default a (random) balanced
    cross-validation is used}
  \item{threshold}{Optional vector of threshold values, for feature selection}
 \item{n.threshold}{Number of threshold value for feature selection (default 10)}



}
\details{
  \code{marginTree.cv} Does K-fold cross-validation of a 
 margin tree classifier, to choose the threshold parameter for
feature selection.}

  
\value{
  A list with components:
\item{threshold}{Threshold values tried}
\item{error}{CV misclassification errors for each threshold}
  \item{size}{Average number of features used at each threshold}
  \item{size.indiv}{Matrix of  number of features used at each split (rows)
 and  at each threshold (columns).}
  \item{yhat}{Matrix of cross-validated predictions.}
  \item{y}{Vector of outcome values}
  \item{folds}{Indices of CV folds}
  \item{call}{The calling sequence}
}

\references{
 Rob Tibshirani and Trevor Hastie. Tech report. Feb. 2006.
Margin trees for high-dimensional classification
  Available at http://www-stat.stanford.edu/~tibs/research.html} 


\author{ Robert Tibshirani and Trevor Hastie}

\examples{
#generate some data with 5 classes and 100 features
set.seed(543)

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

# examine results

cv.obj

}

