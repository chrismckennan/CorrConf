#######This estimates K uses cross validation#######
#It does estimation and testing on the folds in parallel

require(parallel)
require(irlba)

ChooseK_parallel.multB <- function(Y, X=NULL, maxK, B, nFolds=10, A.lin=NULL, c.lin=NULL, D.ker=NULL, Var.0=NULL, tol.rho=1e-3, max.iter.rho=10, svd.method="fast", plotit=T, n_cores=NULL) {
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
  if (is.null(n_cores)) {n_cores <- max(detectCores() - 1, 1)}
  cl <- makeCluster(n_cores)
  clusterEvalQ(cl=cl, {library(irlba); library(CorrConf)})
  clusterExport(cl, c("Y", "B", "maxK", "A.lin", "c.lin", "D.ker", "Var.0", "tol.rho", "max.iter.rho", "svd.method", "folds.rows"), envir=environment())
  out.parallel <- parSapply(cl=cl, 1:nFolds, XVal_K.multB)
  stopCluster(cl)
  
  out$LOO.XV <- 1/n/p * rowSums(out.parallel)
  out$K.hat <- out$K[which.min(out$LOO.XV)]
  if (plotit) {
    plot(out$K, out$LOO.XV, xlab="K", ylab="LOO-XV", main="Leave one out cross validation"); lines(out$K, out$LOO.XV)
    points(out$K.hat, min(out$LOO.XV), pch="x", col="red")
  }

  return(out)
}


######The function to run in parallel######

XVal_K.multB <- function(i) {

  ##Run code##
  Y.1 <- as.matrix(Y[folds.rows != i,]);    #Train with this data (estimate C)
  Y.0 <- as.matrix(Y[folds.rows == i,]);    #Test with these data
  
  train.i <- Optimize.Theta.multB(Y=Y.1, maxK=maxK, B=B, Cov=NULL, A=A.lin, c=c.lin, D.ker=D.ker, Var.0=Var.0, tol.rho=tol.rho, max.iter.rho=max.iter.rho, svd.method=svd.method)
  test.loo.i <- Test.LOOXV.multB(Y.0=Y.0, B=B, train=train.i, A=A.lin, c=c.lin, D.ker=D.ker, Var.0=Var.0)
  
  return(test.loo.i$Loss)
}


####This function performs leave-one-out cross validation on the held out genes after we have estimated C
Test.LOOXV.multB <- function(Y.0, B, train, A=NULL, c=NULL, D.ker=NULL, Var.0=NULL) {
  K <- train$K
  C.list <- train$C
  n <- ncol(Y.0)
  p.0 <- nrow(Y.0)
  SYY.0 <- 1/p.0 * t(Y.0) %*% Y.0
  out <- list()
  out$Loss <- rep(0, length(K))
  
  for (k in K) {
    V.k <- CreateV(B = B, Rho = train$Rho[k+1,])
    V.k <- exp(-1/n * sum(log(svd.wrapper(V.k)$d))) * V.k; sqrt.Vinv.k <- sqrt.mat2(V.k)$Rinv
    if (k == 0) {
      Y.k <- Y.0 %*% sqrt.Vinv.k
      out$Loss[k+1] <- sum(Y.k^2)
    } else {
      C.k <- C.list[[k+1]]
      Y.k <- Y.0 %*% sqrt.Vinv.k; C.k <- sqrt.Vinv.k %*% C.k
      
      SCC.k <- solve(t(C.k) %*% C.k)
      L.k <- Y.k %*% C.k %*% SCC.k
      H.k <- apply(C.k, 1, function(x, Ak) {sum(x * (Ak %*% x))}, Ak=SCC.k)
      R.k <- Y.k - L.k %*% t(C.k)
      out$Loss[k+1] <- sum( sweep(x = R.k, MARGIN = 2, 1/(1-H.k), FUN = "*")^2 )
    }
  }
  return(out)
}
