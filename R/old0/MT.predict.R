MT.predict=function(a,xx,threshold=0){

nbetas<<-NULL
# x is p by n
nsplits=nrow(a$marg.tree$merge)

predfunc=function(a,x1, ii){
  aa=1.5+.5*mypredict(a$svm.splitters[[ii]],t(x1),threshold);
  
   if(a$marg.tree$merge[ii,aa]<0){ pred= -a$marg.tree$merge[ii,aa]}
   else{pred=predfunc(a,x1,a$marg.tree$merge[ii,aa])}
return(pred)
}


n=ncol(xx)

yhat=rep(NA,n)
for(i in 1:n){
#cat(i,fill=T)
  yhat[i]=a$ynams[predfunc(a,xx[,i,drop=F],nsplits)]
}

return(yhat)
}
