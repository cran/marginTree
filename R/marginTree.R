marginTree=function(x,y, method="complete",   n.threshold=20, predict.trainingset=TRUE){

require(e1071)
y=as.character(y)
n=ncol(x)
ynams=unique(y)

# the input x is n by p; make it p by n for convenience
x=t(x)
xs=transf(x)

cat("Computing pairwise margins",fill=T)
marg.obj=compute.marg.matrix(xs,y)

marg.tree=hclust(as.dist(marg.obj$marg),method=method)

junk=compute.splitters(marg.tree,x,xs,y)


return.list=list(marg.obj=marg.obj,marg.tree=marg.tree,svm.splitters=junk$svm.splitters, plot.heights=junk$plot.heights, nclasses=length(table(y)), nlist=junk$nlist, ynams=ynams)

minrat=getmin.ratio(return.list)
#minrat=max(minrat,0)

threshold=NULL
err=NULL
yhat=NULL
nonzero=NULL
if(predict.trainingset){
cat("Computing thresholded estimates",fill=T)
threshold=seq(1,minrat,length=n.threshold)
err=rep(NA,n.threshold)
nonzero=rep(NA,n.threshold)
yhat=matrix(NA,nrow=length(y),ncol=n.threshold)
nonzero=rep(0,n.threshold)
for(i in 1:n.threshold){
cat(c(i,".- Threshold=",round(threshold[i],4)),fill=T)
 yhat[,i]= marginTree.predict(return.list,t(x),threshold[i])
 err[i]=sum(y!=yhat[,i])

 junk=marginTree.getnonzero(return.list, threshold[i])$feature.scores
 for(ii in 1:length(junk)){
   nonzero[i]=nonzero[i]+nrow(junk[[ii]])/length(junk)
  }
}
}


return.list$threshold=threshold
return.list$err=err
return.list$y=y
return.list$yhat=yhat
return.list$nonzero=nonzero
class(return.list)="marginTree"
return(return.list)

}
