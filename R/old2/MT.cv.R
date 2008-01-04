MT.cv <- function(x, y, object, nfold = min(table(y)), folds = NULL, threshold = NULL, n.threshold=10)

{

y=as.character(y)

        this.call <- match.call()
        n <- length(y)

if(is.null(nfold)) {nfold <- min(table(y))}

        if(is.null(folds)) {
                folds <- balanced.folds(y)
        }
         
nfold<- length(folds)

if(missing(threshold)) {
                if(missing(object))
               stop("Must either supply threshold argument, or an MT object")
        }

       
minrat=getmin.ratio(object)
minrat=max(minrat,0)

threshold=seq(1,minrat,length=n.threshold)

        n.threshold <- length(threshold)        ### Set up the data structures
        yhat <- matrix(NA,nrow=n,ncol= n.threshold)
        dimnames(yhat)[[2]] <- paste(seq(n.threshold))
        n.class <- table(y)
        size <- double(n.threshold)

        for(ii in 1:nfold) {
                cat("Fold", ii, ":")
b2=MT.train(x[,-folds[[ii]]],y[-folds[[ii]]],method="complete", predict.trainingset=F)


          for(j in 1:n.threshold){
        yhat[folds[[ii]],j]=MT.predict(b2,x[, folds[[ii]],drop=F],threshold=threshold[j])
        }}

        error <- rep(NA, n.threshold)
        
        for(i in 1:n.threshold) {
               error[i] <- sum(yhat[, i] != y)/n
                }
# find # genes used at each tree split

size.indiv=matrix(NA,nrow=length(object$svm.splitters),ncol=n.threshold)
for(j in 1:n.threshold){
for( i in 1:length(object$svm.splitters)){
  a=object$svm.splitters[[i]]
 rat=a$margshrink[,2]/max(a$margshrink[,2])
rat[1]=1
 o=(1:length(rat))[rat>=threshold[j]]
         o=o[length(o)]
         o[is.na(o)]=length(rat)
         size.indiv[i,j]=a$margshrink[o,1]
}}
size=apply(size.indiv,2,mean)




obj<- list(threshold=threshold, error=error,  size=size, size.indiv=size.indiv, yhat=yhat,y=y,folds=folds, 
                call = this.call)
        class(obj) <- "MTcv"
        obj
}

