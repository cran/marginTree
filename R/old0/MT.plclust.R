MT.plclust=function(train.obj, ...){

marg.tree=train.obj$marg.tree
marg.tree$height=train.obj$plot.heights

plclust(marg.tree,sub="",xlab="",frame.plot=T)

}

