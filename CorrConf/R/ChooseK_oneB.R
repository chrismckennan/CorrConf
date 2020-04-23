##This performs parallelized cross-validation and allows the delta's to be different
#It assumes there is only 1 block matrix B, which is assumed to be a matrix
require(parallel)
require(irlba)

ChooseK_parallel <- function(Y, X=NULL, maxK, B, nFolds=10, tol.rho=1e-3, max.iter.rho=15, svd.method="fast", plotit=T, n_cores=NULL) {
  if (maxK < 1) {
    return(0)
  }

  out <- list()
  out$K <- 0:maxK
  
  if (!is.null(X)) {
    Q <- qr.Q(qr(X), complete=T)[,(ncol(X)+1):nrow(X)]
    B <- t(Q) %*% B %*% Q
    Y <- Y %*% Q
  }
  s.B <- svd.wrapper(B)
  Y <- Y %*% s.B$u
  Lambda <- s.B$d
  
  p <- nrow(Y)
  n <- ncol(Y)
  nFolds <- min(nFolds, n, p)
  
  ##Permute rows of Y##
  perm.rows <- order(runif(p))
  Y <- Y[perm.rows,]
  
  ##Perform K-fold cross validation##
  folds.rows <- cut(1:p, breaks=nFolds, labels=FALSE)
  
  if (is.null(n_cores)) {n_cores <- max(detectCores() - 1, 1)}
  cl <- makeCluster(n_cores)
  clusterEvalQ(cl=cl, {library(irlba); library(CorrConf)})
  clusterExport(cl, c("Y", "Lambda", "maxK", "folds.rows", "tol.rho", "max.iter.rho", "svd.method"), envir=environment())
  out.parallel <- parSapply(cl=cl, 1:nFolds, XVal_K)
  stopCluster(cl)
  
  out$LOO.XV <- 1/n/p * apply(out.parallel, 1, sum)
  out$K.hat <- out$K[which.min(out$LOO.XV)]
  if (plotit) {
    plot(out$K, out$LOO.XV, xlab="K", ylab="LOO-XV", main="Leave one out cross validation"); lines(out$K, out$LOO.XV)
    points(out$K.hat, min(out$LOO.XV), pch="x", col="red")
  }
  
  return(out)
}




######Routine within cross-validation loop######
#i is current fold, which is the variable that will be parallelized. All other variables need to be passed to the function
XVal_K <- function(i) {
  
  ##Run code##
  Y.1 <- Y[folds.rows != i,];    #Train with this data (estimate C)
  Y.0 <- Y[folds.rows == i,];    #Test with these data
  
  train.i <- Optimize.Theta(Y=Y.1, K=maxK, B=Lambda, Cov=NULL, tol.rho=tol.rho, max.iter.rho=max.iter.rho, svd.method=svd.method, return.C=T)
  test.loo.i <- Test.LOOXV.simrho(Y.0=Y.0, Lambda=Lambda, train=train.i)  ##Test with all delta's being the same, since this tends to minimize the loss function
  
  return(test.loo.i$Loss)
}



##########Subroutines to be used in the above functions##########

##Leave one out cross validated error
Test.LOOXV <- function(Y.0, Lambda, train) {
  K <- train$K
  C.list <- train$C
  n <- ncol(Y.0)
  p.0 <- nrow(Y.0)
  out <- list()
  out$Loss <- rep(0, length(K))
  
  for (k in K) {
    if (k == 0) {
      rho.k <- optimize.rho.K0(Y = Y.0, eigs.K = Lambda, rho.start = seq(0, 0.95, by=0.05))$rho
      V.k <- 1 + rho.k * (Lambda - 1)
      out$Loss[k+1] <- sum(sweep(Y.0, MARGIN = 2, 1/sqrt(V.k), "*")^2)
      
      rho.previous <- rho.k; rho.previous <- max(rho.previous, 0.05)
    } else {
      C.k <- C.list[[k+1]]
      Q.k <- qr.Q(qr(C.k), complete=T)[,(k+1):n]
      s.k <- svd.wrapper(t(Q.k * Lambda) %*% Q.k)
      rho.k <- optimize.rho.K0(Y = Y.0 %*% Q.k %*% s.k$u, eigs.K = s.k$d, rho.start = rho.previous)$rho
      rho.previous <- rho.k; rho.previous <- max(rho.previous, 0.05)
      
      V.k <- 1 + rho.k * (Lambda - 1)
      Y.k <- sweep(Y.0, MARGIN = 2, 1/sqrt(V.k), "*")
      C.k <- C.k / sqrt(V.k)
      L.k <- Y.k %*% C.k %*% solve(t(C.k) %*% C.k)
      H.k <- apply(C.k, 1, function(x, A) {sum(x * (A %*% x))}, A=solve(t(C.k) %*% C.k))
      R.k <- Y.k - L.k %*% t(C.k)
      out$Loss[k+1] <- sum( sweep(x = R.k, MARGIN = 2, 1/(1-H.k), FUN = "*")^2 )
    }
  }
  return(out)
}


