MT.confusion <- function(fit, threshold, extra = TRUE) {
  temp=abs(fit$threshold -threshold)
  ii <- (1:length(fit$threshold))[temp==min(temp)]
  ii <- ii[1]
  predicted <- fit$yhat[, ii]
  
  true <- fit$y
  tt <- table(true, predicted)

  if (extra) {
    tt1 <- tt
    diag(tt1) <- 0
    tt <- cbind(tt, apply(tt1, 1, sum)/apply(tt, 1, sum))
    dimnames(tt)[[2]][ncol(tt)] <- "Class Error rate"
    print(tt)
    cat(c("Overall error rate=", round(sum(tt1)/sum(tt), 3)),
        fill= TRUE)
  }
  if (!extra) {
    return(tt)
  }
}

which.is.max <- function(x)
{
        y <- seq(length(x))[x == max(x)]
        if(length(y) > 1)
                y <- sample(y, 1)
        y
}

