MT.listfeatures=function(train.obj, threshold){

marg.tree=train.obj$marg.tree
marg.tree$height=train.obj$plot.heights

plclust(marg.tree,sub="",xlab="",frame.plot=T)

}

