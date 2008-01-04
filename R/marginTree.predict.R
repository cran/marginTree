marginTree.predict=function(train.obj,x,threshold=1){

require(e1071)
nbetas<<-NULL
# x is n by p; make it p by n for convenience
x=t(x)

nsplits=nrow(train.obj$marg.tree$merge)

predfunc=function(a,x1, ii){
  aa=1.5+.5*mypredict(train.obj$svm.splitters[[ii]],t(x1),threshold);
  
   if(train.obj$marg.tree$merge[ii,aa]<0){ pred= -train.obj$marg.tree$merge[ii,aa]}
   else{pred=predfunc(train.obj,x1,train.obj$marg.tree$merge[ii,aa])}
return(pred)
}


n=ncol(x)

yhat=rep(NA,n)
for(i in 1:n){
  yhat[i]=train.obj$ynams[predfunc(a,x[,i,drop=F],nsplits)]
}

return(yhat)
}
