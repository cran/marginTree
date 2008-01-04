MT.get.nonzero=function(train.obj,threshold){
nn=nrow(train.obj$marg.tree$merge)
feature.scores=vector("list",nn)
left.features=vector("list",nn)
right.features=vector("list",nn)

for(i in 1:nn){
 a=train.obj$svm.splitters[[i]]
        bb=a$beta
        rat = a$margshrink[, 2]/max(a$margshrink[, 2])
        o = (1:length(rat))[rat >= threshold]
        o = o[length(o)]
        o[is.na(o)] = length(rat)
        temp = a$margshrink[o, 1]
        r = rank(-abs(bb))
        bb[r > temp] = 0 
        bb = bb/sqrt(sum(bb^2))
        ord=order(-abs(bb[bb!=0]))
feature.scores[[i]]=cbind((1:length(bb))[bb!=0], bb[bb!=0])[ord,]
dimnames(feature.scores[[i]])=list(NULL,c("Feature#", "Weight"))

left.features[[i]]=my.descendants(train.obj$marg.tree$merge,i,1)
right.features[[i]]=my.descendants(train.obj$marg.tree$merge,i,2)
    }

return.list=list(feature.scores=feature.scores,left.features=left.features,
right.features=right.features)
class(return.list)="MTgetnonzero"
return(return.list)
}

