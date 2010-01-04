\name{marginTree.getnonzero}
\alias{marginTree.getnonzero}
\title{Get important features at each split of the margin tree}
\description{
A function to get important features at each split of the margin tree
}
\usage{
marginTree.getnonzero(train.obj,  threshold)
}

\arguments{
\item{train.obj}{Output from call to marginTree}
\item{threshold}{Threshold for feature selection: between 0 and 1}
}
\details{
  \code{marginTree.getnonzero} Does hard thresholding of the
weight vector at each split in the margin tree, to select features}
\value{
A list-- one element per split in the tree--  with components
\item{feature.scores}{The coefficient for the  selected feature}
\item{left.classes}{The outcome classes assigned to the left branch}
\item{right.classes}{The outcome classes assigned to the right branch}
}

\references{
 Rob Tibshirani and Trevor Hastie. Tech report. Feb. 2006.
Margin trees for high-dimensional classification
  Available at http://www-stat.stanford.edu/~tibs/research.html} 


\author{Rob Tibshirani and Trevor Hastie}

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


# examine the selected features at each split

junk<- marginTree.getnonzero(train.obj,threshold=.5)

summary(junk)
}

