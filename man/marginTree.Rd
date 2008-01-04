\name{marginTree}
\alias{marginTree}
\title{ A function to train  a
 margin tree  classifier}
\description{
A function to train a
 margin tree  classifier. This is a hierachical version of the support
vector classifier, useful for more than 2 classes.
}
\usage{
marginTree(x,y, method="complete",  n.threshold=20, predict.trainingset=TRUE)
}

\arguments{
  \item{x}{The input data of feature values, n samples by p features}
  \item{y}{Class labels- vector of length n}
  \item{method}{Clustering method- "complete" (default and recommended);  "average" or "single"}
  \item{n.threshold}{Number of threshold values desired (default 10)}
\item{predict.trainingset}{Predict the training set?
  Required for computing error rates and in preparation for cross-validation.
Slows the computation down a bit for large datasets.
Default TRUE}
}
\details{
  \code{marginTree} fits a margin tree classifier.
It is useful with more than 2 outcome classes, and when the
number of features exceeds the number of observations, for example
in genomic and proteomic applications.
  Details may be found in the paper on the website listed below.
This function calls the \code{svm} function from library \code{e1071}.
}

\value{
  A list with components:
  \item{marg.obj}{List with components \code{marg}- the matrix of pairwise
margins, and \code{svmfit}, a list of svm classifiers used to 
find each pairwise margin}
  \item{marg.tree}{The margin tree. Same format as that produced by  hclust.}
  \item{svm.splitters}{List a svm objects for classification at each
 junction in the margin tree.}
  \item{plot.heights}{Heights for plotting the margin tree}
  \item{nclasses}{Number of outcome classes}
  \item{nlist}{For internal use}
\item{ynams}{Character names of outcome classes}
\item{threshold}{Threshold values tried}
\item{err}{Number of training errors for each threshold value}
\item{y}{Training class labels}
\item{yhat}{Matrix of predicted classes, one column per threshold value}
\item{nonzero}{Average number of  nonzero features per threshold value}
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



#train the classifier

train.obj<-   marginTree(x,y)

# print out results to see training error rates
 
  train.obj
}
\keyword{ }

