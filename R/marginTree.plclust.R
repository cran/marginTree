marginTree.plclust=function(train.obj, ...){
if(train.obj$nclasses==2){
stop("No cluster tree is  available since there are only 2 classes in the data")
}
else{
marg.tree=train.obj$marg.tree
marg.tree$height=train.obj$plot.heights

plclust(marg.tree,sub="",xlab="",frame.plot=T)
}

}

