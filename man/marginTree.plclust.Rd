\name{marginTree.plclust}
\alias{marginTree.plclust}
\title{Plot a  margin tree}
\description{
A function to plot margin tree}
}
\usage{
marginTree.plclust(train.obj, ...)
}

\arguments{
\item{train.obj}{Output from call to marginTree}
\item{...}{Additional arguments to pass to the function plclust}
}

\details{
  \code{marginTree.plclust}{Plots a margin tree,
with height of the vertical gap between each tree junction equal
to the achieved margin for that split}.
}

  
\value{
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


# plot the margin tree
marginTree.plclust(train.obj)

}
\keyword{ }

