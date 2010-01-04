\name{marginTree.confusion}
\alias{marginTree.confusion}
\title{ A function giving a table of true versus predicted values,
  from a margin tree.}
\description{A function giving a table of true versus predicted values,
  from a margin tree}
\usage{
marginTree.confusion(train.obj, threshold, extra=TRUE)
}

\arguments{
  
  \item{train.obj}{The result of a call to marginTree.train or marginTree.cv}
  \item{threshold}{The desired threshold value}
\item{extra}{Should the  classwise and overall error rates be returned?
Default TRUE}
  
  
}    

\details{
  \code{marginTree.confusion} Gives a cross-tabulation of true versus
  predicted classes for the train.obj returned by marginTree or marginTree.cv,
  at the specified threshold.
}






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

train.obj<-   marginTree(x,y)

# print out confusion matrix 

marginTree.confusion(train.obj, threshold=.5)

# do cross-validation and print out CV confusion matrix
cv.obj<-marginTree.cv(x,y,train.obj)

marginTree.confusion(cv.obj, threshold=.5)
}

