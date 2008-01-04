marginTree.plotcv <- function(cv.obj) {
  par(mar = c(5, 5, 5, 1))
  par(mfrow = c(2, 1))
  n <- nrow(cv.obj$yhat)
  y <- cv.obj$y
  nc <- length(table(y))
  nfolds <- length(cv.obj$folds)
  err <- matrix(NA, ncol = ncol(cv.obj$yhat), nrow = nfolds)
  temp <- matrix(y, ncol = ncol(cv.obj$yhat), nrow = n)
  ni <- rep(NA, nfolds)
  for(i in 1:nfolds) {
    ii <- cv.obj$folds[[i]]
    ni[i] <- length(cv.obj$folds[[i]])
    err[i,  ] <- apply(temp[ii,  ] != cv.obj$yhat[ii,  ], 2, sum)/ni[i]
  }
  se <- sqrt(apply(err, 2, var)/nfolds)
  plot(cv.obj$threshold, cv.obj$error, ylim = c(-0.1, 0.8), xlab = 
       "Value of threshold  ", ylab = "Misclassification Error", type
       = "n", yaxt = "n")
  axis(3, at = cv.obj$threshold, lab = paste(cv.obj$size), srt = 90, adj = 0)
  mtext("Number of features", 3, 4, cex = 1.2)
  axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8))
  lines(cv.obj$threshold, cv.obj$error, col = 2)
  o <- cv.obj$err == min(cv.obj$err)
  points(cv.obj$threshold[o], cv.obj$error[o], pch = "x")
  error.bars(cv.obj$threshold, cv.obj$err - se, cv.obj$err + se)
  err2 <- matrix(NA, nrow = length(unique(y)), ncol = length(cv.obj$threshold
                                                 ))
  for(i in 1:(length(cv.obj$threshold) )) {
    s <- marginTree.confusion(cv.obj, cv.obj$threshold[i], extra = FALSE)
    diag(s) <- 0
    err2[, i] <- apply(s, 1, sum)/table(y)
  }
  plot(cv.obj$threshold, err2[1,  ],xlim=c(min(cv.obj$threshold),1.3), ylim = c(-0.1, 1.1), xlab = 
       "Value of threshold ", ylab = "Misclassification Error", type
       = "n", yaxt = "n")
  axis(3, at = cv.obj$threshold, lab = paste(cv.obj$size), srt = 90, adj = 0)     
  axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8))
  for(i in 1:nrow(err2)) {
    lines(cv.obj$threshold, err2[i,  ], col = i + 1)
  }
  legend(1.1, 1.1, dimnames(table(y))[[1]], col = (2:(nc + 1)), lty = 1,cex=.6)
  par(mfrow = c(1, 1))
}

error.bars <-function(x, upper, lower, width = 0.02, ...) {
  xlim <- range(x)
  barw <- diff(xlim) * width
  segments(x, upper, x, lower, ...)
  segments(x - barw, upper, x + barw, upper, ...)
  segments(x - barw, lower, x + barw, lower, ...)
  range(upper, lower)
}

