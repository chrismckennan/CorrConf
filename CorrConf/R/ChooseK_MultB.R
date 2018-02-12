#######This estimates K uses cross validation#######
#It does estimation and testing on the folds in parallel

require(parallel)
require(irlba)

ChooseK_parallel.multB <- function(Y, X=NULL, maxK, B, nFolds=10, tol.rho=1e-3, max.iter.rho=10, svd.method="fast") {
  if (maxK < 1) {
    return(0)
  }
  
  out <- list()
  out$K <- 0:maxK
  
  if (!is.null(X)) {
    Q <- qr.Q(qr(X), complete=T)[,(ncol(X)+1):nrow(X)]
    if (!is.null(B)) {B <- lapply(B, function(x, Q){t(Q) %*% x %*% Q}, Q=Q)}
    Y <- Y %*% Q
  }
  
  p <- nrow(Y)
  n <- ncol(Y)
  nFolds <- min(nFolds, n, p)
  
  ##Permute rows of Y##
  perm.rows <- order(runif(p))
  Y <- Y[perm.rows,]
  
  ##Perform K-fold cross validation##
  folds.rows <- cut(1:p, breaks=nFolds, labels=FALSE)
  
  n_cores <- max(detectCores() - 1, 1)
  cl <- makeCluster(n_cores)
  clusterEvalQ(cl=cl, {library(irlba); library(CorrConf)})
  clusterExport(cl, c("Y", "B", "maxK", "tol.rho", "max.iter.rho", "svd.method", "folds.rows"), envir=environment())
  out.parallel <- parLapply(cl=cl, 1:nFolds, XVal_K.multB)
  stopCluster(cl)
  
  out$LOO.XV <- 1/n/p * apply(out.parallel, 1, sum)
  out$K.hat <- out$K[which.min(out$LOO.XV)]
  plot(out$K, out$LOO.XV, xlab="K", ylab="LOO-XV", main="Leave one out cross validation"); lines(out$K, out$LOO.XV)
  points(out$K.hat, min(out$LOO.XV), pch="x", col="red")
  
  return(out)
}


######The function to run in parallel######

XVal_K.multB <- function(i) {

  ##Run code##
  Y.1 <- Y[folds.rows != i,];    #Train with this data (estimate C)
  Y.0 <- Y[folds.rows == i,];    #Test with these data
  
  train.i <- Optimize.Theta.multB(Y=Y.1, maxK=maxK, B=B, Cov=NULL, tol.rho=tol.rho, max.iter.rho=max.iter.rho, svd.method="fast")
  test.loo.i <- Test.LOOXV.multB(Y.0=Y.0, B=B, train=train.i)
  
  return(list(Rho=test.i$Rho, Loss=test.i$Fnorm2, Loss.LOO=test.loo.i$Loss, CorrE=test.loo.i$CorrE, Loss.noCorr=test.nocorr.i$Loss))  
}


####This function performs leave-one-out cross validation on the held out genes after we have estimated C
Test.LOOXV.multB <- function(Y.0, B, train) {
  Rho <- NULL
  K <- train$K
  C.list <- train$C
  n <- ncol(Y.0)
  p.0 <- nrow(Y.0)
  SYY.0 <- 1/p.0 * t(Y.0) %*% Y.0
  out <- list()
  out$Loss <- rep(0, length(K))
  
  for (k in K) {
    if (k == 0) {
      if (is.null(Rho)) {
        out.rho.k <- Est.Corr.multB(Y=SYY.0, B=B, simple.rho=T)
        Rho.k <- out.rho.k$Rho; Rho.previous <- Rho.k; Rho.previous[Rho.previous < 0.05] <- 0.05
      }
      V.k <- CreateV(B = B, Rho = Rho.k); sqrt.Vinv.k <- sqrt.mat2(V.k)$Rinv
      Y.k <- Y.0 %*% sqrt.Vinv.k

      out$Loss[k+1] <- sum(Y.k^2)
    } else {
      C.k <- C.list[[k+1]]
      Q.k <- qr.Q(qr(C.k), complete=T)[,(k+1):n]
      if (is.null(Rho)) {
        out.rho.k <- Est.Corr.multB(Y = t(Q.k) %*% SYY.0 %*% Q.k, B = lapply(B, function(x, Q.k){t(Q.k) %*% x %*% Q.k}, Q=Q.k), theta.0 = Rho.previous, simple.rho=T)
        Rho.k <- out.rho.k$Rho; Rho.previous <- Rho.k; Rho.previous[Rho.previous < 0.05] <- 0.05
      } else {
        Rho.k <- Rho
      }
      V.k <- CreateV(B = B, Rho = Rho.k); sqrt.Vinv.k <- sqrt.mat2(V.k)$Rinv
      Y.k <- Y.0 %*% sqrt.Vinv.k; C.k <- sqrt.Vinv.k %*% C.k

      L.k <- Y.k %*% C.k %*% solve(t(C.k) %*% C.k)
      H.k <- apply(C.k, 1, function(x, A) {sum(x * (A %*% x))}, A=solve(t(C.k) %*% C.k))
      R.k <- Y.k - L.k %*% t(C.k)
      out$Loss[k+1] <- sum( sweep(x = R.k, MARGIN = 2, 1/(1-H.k), FUN = "*")^2 )
    }
  }
  return(out)
}
