marginTree.getnonzero=function(train.obj,threshold){
nn=nrow(train.obj$marg.tree$merge)
feature.scores=vector("list",nn)
left.classes=vector("list",nn)
right.classes=vector("list",nn)

ynams=train.obj$ynams

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

left.classes[[i]]=ynams[my.descendants(train.obj$marg.tree$merge,i,1)]
right.classes[[i]]=ynams[my.descendants(train.obj$marg.tree$merge,i,2)]
    }

return.list=list(feature.scores=feature.scores,left.classes=left.classes,
right.classes=right.classes)
class(return.list)="marginTreegetnonzero"
return(return.list)
}

