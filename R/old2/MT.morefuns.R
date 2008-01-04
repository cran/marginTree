"mysvm"<-function(x,y,K, alpha0=rep(0,len(y)), deg=NULL, cxx=NULL){
#solves  sv problem
# K is C in Vapnik page 133
# returns normal to plane w0 and support coefs alpha

save.call<-match.call()

n<-nrow(x)
p<-ncol(x)
nclin<-1

cc<-y
alpha0<-rep(0,n)

storage.mode(cc)<-"double"
BIG<-10e10

bl<-c(rep(0,n), 0)
bu<-c(rep(K,n),0)
cvec<-rep(-1,n)
y<-matrix(y,ncol=1)

if(is.null(deg) & is.null(cxx)){cxx=x%*%t(x)}
if(is.null(deg)){ a<-.5*(y%*%t(y))*cxx;  a<-(a+t(a))/2}
if(!is.null(deg)){
  a<-ipf(x,x,deg)  
  a<-.5*(y%*%t(y))*a
  a<-(a+t(a))/2
}
# not used

istate<-rep(0, n+nclin)
kx<-rep(0,n)

leniw<-n
lenw<-2*n*n+10*n+6
storage.mode(cc)<-"double"
storage.mode(bl)<-"double"
storage.mode(bu)<-"double"
storage.mode(alpha0)<-"double"
storage.mode(x)<-"double"
storage.mode(y)<-"double"

storage.mode(cvec)<-"double"
storage.mode(a)<-"double"
storage.mode(istate)<-"integer"
storage.mode(leniw)<-"integer"
storage.mode(lenw)<-"integer"
dyn.load("/home/tibs/PAPERS/sv/lssol/LSSOL.so")

junk<-.Fortran("mysolve",
 as.integer(n),
 as.integer(n),
 as.integer(nclin),
 as.integer(nclin),
 as.integer(n),
 as.matrix(cc), 
 as.double(bl),
 as.double(bu),
 as.double(cvec),
 istate=as.integer(istate),
 as.integer(kx), 
 alpha=as.double(alpha0),
 as.matrix(a),
 double(n),
 inform=integer(1),
 iter=integer(1), 
 obj=double(1),
 clamda=double(n+1),
 iw=integer(leniw),
 as.integer(leniw),
 w=double(lenw), 
as.integer(lenw))

w0<-t(x)%*%(junk$alpha*y)


yhat0=x%*%w0/2
if(mean(yhat0[y==1]) > mean(yhat0[y==-1])) {
   gap=min(yhat0[y==1])-max(yhat0[y==-1])
   cutp=max(yhat0[y==-1])+gap/2
   dir=+1
}
else{
   gap=min(yhat0[y==-1])-max(yhat0[y==1])
   cutp=max(yhat0[y==1])+gap/2
   dir=-1
}

return(list(beta=w0/2, coefs=junk$alpha[junk$alpha>0],crit=junk$obj, info=junk$inform, iter=junk$iter, istate=junk$istate,supp=sum(junk$alpha>0)/(length(junk$alpha)-1), deg=deg, y=y,gap=gap,dir=dir,cutp=cutp,index=(1:length(y))[junk$alpha>0],
call=save.call))
}


mypredict=function(a,x,threshold){
# prediction with thresholded betas
# takes object "a" from compute.splitters
beta=a$beta

vars=1:ncol(x)
meanx=rep(0,ncol(x))
mu=0
cutp=a$cutp
dir=a$dir

if(threshold>0){
 rat=a$margshrink[,2]/max(a$margshrink[,2])
         o=(1:length(rat))[rat>=threshold]
         o=o[length(o)]
         o[is.na(o)]=length(rat)
         temp=a$margshrink[o,1]
         r=rank(-abs(beta))
         beta[r>temp]=0
         beta=beta/sqrt(sum(beta^2))
         cutp=a$cutp.shrink[o]
         dir=a$dir.shrink[o]
}


nbetas<<-c(nbetas,sum(beta!=0))
yhat=scale(x[,vars,drop=F],meanx,FALSE)%*%beta +mu
if(dir==1){ yhatt=-1+2*(yhat>cutp)}
if(dir== -1){yhatt=1-2*(yhat>cutp)}
return(yhatt)
}
margfunc=function(x,cxx,y,a=NULL){
if(is.null(a)){
  a=mysvm(t(x),y,K=10e99, cxx=cxx)
#  a=svm(t(x),as.factor(y),cost=10e99, scale=F,  kernel="linear")
}
beta=a$beta
#beta = x[, a$index] %*% a$coefs
margin=2/sqrt(sum(beta^2))
return(list(margin=margin, svmfit=a))
}


getmin.ratio=function(a){
nn=length(a$svm.splitters)
rat=NULL
for(i in 1:nn){
val=a$svm.splitters[[i]]$margshrink[,2]
rat=c(rat,min(val)/max(val))
}
return(min(rat))
}


compute.splitters=function(marg.tree,x,cxx,y, compute.ngenes=TRUE ){
cat("Computing splitters",fill=T)
n=length(y)
nams=unique(y)
nsplits=nrow(marg.tree$merge)
nlist=NULL
svm.splitters=vector("list",nsplits)
for(i in 1:nsplits){
cat(c("Split=",i),fill=T)
   o1=my.descendants(marg.tree$merge,i,1)
   o2=my.descendants(marg.tree$merge,i,2)
   oo1=NULL;oo2=NULL
   for(ii  in 1:length(o1)){
     oo1=c(oo1,(1:n)[y==nams[o1[ii]]])
   }
   for(ii in 1:length(o2)){
     oo2=c(oo2,(1:n)[y==nams[o2[ii]]])
   }

   xx=x[,c(oo1,oo2)]
   yy=c(rep(-1,length(oo1)),rep(1,length(oo2)))
   cxxtemp=cxx[c(oo1,oo2),c(oo1,oo2)]
   svm.splitters[[i]]=mysvm(t(xx),yy,K=10e99,cxx=cxxtemp)
#svm.splitters[[i]]=svm(t(xx),as.factor(yy),cost=10e99, scale=F,  kernel="linear")

marg=margfunc(xx,cxxtemp,as.factor(yy),svm.splitters[[i]])$marg
a=svm.splitters[[i]]
#beta=xx[,a$index]%*%a$coefs
beta=a$beta

svm.splitters[[i]]$left.classes=o1
svm.splitters[[i]]$right.classes=o2
svm.splitters[[i]]$left.ind=oo1
svm.splitters[[i]]$right.ind=oo2

if(compute.ngenes){

 r=rank(-abs(beta))
nlist=trunc(exp(seq(log(nrow(x)),log(2),length=50)))
nlist=nlist[nlist>10]

margshrink=matrix(NA,nrow=length(nlist),ncol=2)

ii=0
dir.shrink=NULL
cutp.shrink=NULL
for(ngenes in nlist){
temp=beta
oo=r<=ngenes
temp=temp[oo]
temp=temp/sqrt(sum(temp^2))
yhat=t(xx[oo,])%*%temp
ii=ii+1
if(mean(yhat[yy==1]) > mean(yhat[yy==-1])) {
   gap=min(yhat[yy==1])-max(yhat[yy==-1])
   cutp.shrink=c(cutp.shrink,max(yhat[yy==-1])+gap/2)
   dir.shrink=c(dir.shrink,1)
}
else{
   gap=min(yhat[yy==-1])-max(yhat[yy==1])
   cutp.shrink=c(cutp.shrink,max(yhat[yy==1])+gap/2)
   dir.shrink=c(dir.shrink,-1)
}
margshrink[ii,]=c(ngenes,gap)


}
svm.splitters[[i]]$margshrink=margshrink
svm.splitters[[i]]$cutp.shrink=cutp.shrink
svm.splitters[[i]]$dir.shrink=dir.shrink
}


yhat0=t(xx)%*%beta

if(mean(yhat0[yy==1]) > mean(yhat0[yy==-1])) {
   gap=min(yhat0[yy==1])-max(yhat0[yy==-1])
   cutp=max(yhat0[yy==-1])+gap/2
   dir=+1
}
else{
   gap=min(yhat0[yy==-1])-max(yhat0[yy==1])
   cutp=max(yhat0[yy==1])+gap/2
   dir=-1
 }
 svm.splitters[[i]]$marg=marg
 svm.splitters[[i]]$cutp=cutp
 svm.splitters[[i]]$dir=dir
svm.splitters[[i]]$beta=beta
}
# compute heights for plotting tree
marg=NULL
for(i in 1:length(svm.splitters)){
marg=c(marg,svm.splitters[[i]]$marg)
}

plot.heights=rep(0,length(marg))
for(i  in length(marg):1){
  o1=marg.tree$merge[i,1]
  o2=marg.tree$merge[i,2]
if(o1>0){plot.heights[o1]=plot.heights[i]-marg[i]}
if(o2>0){plot.heights[o2]=plot.heights[i]-marg[i]}

}

plot.heights=plot.heights-min(plot.heights)

return(list(marg.tree=marg.tree,svm.splitters=svm.splitters, plot.heights=plot.heights, nclasses=length(table(y)), nlist=nlist))
}



 descendants <- function(m,k){
  # the done object indicates what rows of m were used
    done <- k
    if (m[k,1]<0) left <- -m[k,1]
    else {
      junk <- descendants(m,m[k,1])
      left <- junk[[1]]
      done <- c(done,junk[[2]])
    }
    if (m[k,2]<0) right <- -m[k,2]
    else {
      junk <- descendants(m,m[k,2])
      right <- junk[[1]]
      done <- c(done,junk[[2]])
    }
    return(list(c(left,right),done))
  }
my.descendants<- function(m,j,k){
  #returns descendants of merge element m[j,k]
  #uses a kluge in order to make us of "descendants" function

  if(k==1){m[j,2]<- -99999}
  if(k==2){m[j,1]<- -99999}
junk<- sort(descendants(m,j)[[1]])
return(junk[-length(junk)])
}
compute.marg.matrix=function(x,cxx,y){
k=length(table(y))

n=ncol(x)

svmfit=vector("list",k*(k-1)/2)
marg=matrix(NA,nrow=k,ncol=k)
nams=unique(y)
ii=0
for(i in 1:(k-1)){
for(j in (i+1):k){
cat(c(i,j),fill=T)
o=y==nams[i] | y==nams[j]
xx=x[,o]
cxxtemp=cxx[o,o]
yy= 1-2*(y[o]==nams[i])
junk=margfunc(xx,cxxtemp,yy)
marg[i,j]=junk$margin
ii=ii+1
svmfit[[ii]]=junk
dimnames(marg)=list(nams,nams)
}}
marg[is.na(marg)]=0
marg=marg+t(marg)
return(list(marg=marg,svmfit=svmfit))
}

balanced.folds <- function(y, nfolds = min(min(table(y)), 10)) {
   totals <- table(y)
   fmax <- max(totals)
   nfolds <- min(nfolds, fmax)
   nfolds= max(nfolds, 2)
         # makes no sense to have more folds than the max class size
   folds <- as.list(seq(nfolds))
   yids <- split(seq(y), y)
         # nice we to get the ids in a list, split by class
###Make a big matrix, with enough rows to get in all the folds per class
   bigmat <- matrix(NA, ceiling(fmax/nfolds) * nfolds, length(totals))
   for(i in seq(totals)) {
     if(length(yids[[i]])>1){bigmat[seq(totals[i]), i] <- sample(yids[[i]])}
     if(length(yids[[i]])==1){bigmat[seq(totals[i]), i] <- yids[[i]]}

   }
   smallmat <- matrix(bigmat, nrow = nfolds)# reshape the matrix
### Now do a clever sort to mix up the NAs
   smallmat <- permute.rows(t(smallmat))   ### Now a clever unlisting
         # the "clever" unlist doesn't work when there are no NAs
         #       apply(smallmat, 2, function(x)
         #        x[!is.na(x)])
   res <-vector("list", nfolds)
   for(j in 1:nfolds) {
     jj <- !is.na(smallmat[, j])
     res[[j]] <- smallmat[jj, j]
   }
   return(res)
 }
permute.rows <-function(x)
{
        dd <- dim(x)
        n <- dd[1]
        p <- dd[2]
        mm <- runif(length(x)) + rep(seq(n) * 10, rep(p, n))
        matrix(t(x)[order(mm)], n, p, byrow = TRUE)
}

print.MTtrain <- function(x, ...) {
  cat("Call:\n")
  dput(x$call)
  mat <- rbind(threshold = format(round(x$threshold, 3)), nonzero = 
               format(trunc(x$nonzero)), errors = x$err)
  dimnames(mat) <- list(dimnames(mat)[[1]], paste(1:ncol(mat)))
  print(t(mat), quote = FALSE)
  invisible()
}


 print.MTcv <-function(x, ...) {
   cat("Call:\n")
   dput(x$call)
  
   mat <- rbind(threshold = format(round(x$threshold, 3)), nonzero = 
                format(trunc(x$size)), errors = trunc(x$err * nrow(
                                              x$yhat)))
   dimnames(mat) <- list(dimnames(mat)[[1]], paste(1:ncol(mat)))
 
 cat("\n")
   print(t(mat), quote = FALSE)
   invisible()
}

