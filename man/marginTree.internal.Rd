\name{marginTree-internal}
\title{Internal marginTree functions}
\alias{mypredict}
\alias{mysvd}
\alias{margfunc}
\alias{getmin.ratio}
\alias{compute.splitters}
\alias{balanced.folds}
\alias{permute.rows}
\alias{print.marginTree}
\alias{print.marginTreecv}
\alias{print.marginTreegetnonzero}
\alias{summary.marginTreegetnonzero}
\alias{transf}
\description{Internal marginTree functions}
\usage{
mypredict(a,x,threshold)
mysvd(x, n.components = NULL) 
margfunc(x,y,a=NULL)
getmin.ratio(a)
compute.splitters(marg.tree,x,xs,y, compute.ngenes=TRUE )
balanced.folds(y, nfolds = min(min(table(y)), 10))
permute.rows(x)
print.marginTree(x, ...) 
print.marginTreecv(x, ...) 
print.marginTreegetnonzero(x, ...) 
summary.marginTreegetnonzero(object, ...) 
transf(x,xx=NULL)
}
\author{Rob Tibshirani and Trevor Hastie}
\details{
  These are not to be called by the user.
}
\keyword{internal}


