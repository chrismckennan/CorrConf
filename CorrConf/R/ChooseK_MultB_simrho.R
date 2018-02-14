#######This estimates K uses cross validation#######
#It does estimation and testing on the folds in parallel

require(parallel)
require(irlba)

ChooseK_parallel.multB.simrho <- function(Y, X=NULL, maxK, B, nFolds=10, tol.rho=1e-3, max.iter.rho=10, svd.method="fast", plotit=T) {
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
  
  SYY <- 1/p * t(Y) %*% Y
  
  ##Permute rows of Y##
  perm.rows <- order(runif(p))
  Y <- Y[perm.rows,]
  
  ##Perform K-fold cross validation##
  folds.rows <- cut(1:p, breaks=nFolds, labels=FALSE)
  Y.list <- vector(mode = "list", length = nFolds)    #Folds of Y as a list
  for (f in 1:nFolds) {
    Y.list[[f]] <- Y[folds.rows==f,]
  }
  rm(Y)
  
  n_cores <- max(detectCores() - 1, 1)
  cl <- makeCluster(n_cores)
  clusterEvalQ(cl=cl, {library(irlba); library(CorrConf)})
  clusterExport(cl, c("SYY", "B", "maxK", "tol.rho", "max.iter.rho", "svd.method", "folds.rows", "p"), envir=environment())
  out.parallel <- parSapply(cl=cl, Y.list, XVal_K.multB.simrho)
  stopCluster(cl)
  
  out$LOO.XV <- 1/n/p * apply(out.parallel, 1, sum)
  out$K.hat <- out$K[which.min(out$LOO.XV)]
  if (plotit) {
    plot(out$K, out$LOO.XV, xlab="K", ylab="LOO-XV", main="Leave one out cross validation"); lines(out$K, out$LOO.XV)
    points(out$K.hat, min(out$LOO.XV), pch="x", col="red")
  }
  
  return(out)
}


######The function to run in parallel######

XVal_K.multB.simrho <- function(Y.test) {
  p.1 <- nrow(Y.test)
  p.0 <- p - p.1
  
  SYY.0 <- 1/p.0 * (p * SYY - t(Y.test) %*% Y.test)
  train.i <- Optimize.Theta.multB.simrho(SYY = SYY.0, B = B, maxK = maxK, tol.rho = tol.rho, max.iter.rho = max.iter.rho, svd.method = svd.method)
  test.loo.i <- Test.LOOXV.multB(Y.0=Y.test, B=B, train=train.i)
  return(test.loo.i$Loss)
}
