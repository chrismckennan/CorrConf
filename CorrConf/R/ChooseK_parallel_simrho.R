require(parallel)

##This chooses K

ChooseK_parallel.simrho <- function(Y, X=NULL, maxK, B, nFolds, tol.rho=1e-3, max.iter.rho=15, svd.method="fast", plotit=T, n_cores=NULL) {
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
  SYY <- 1/p * t(Y) %*% Y
  
  ##Permute rows of Y##
  Y <- Y[order(runif(p)),]
  
  ##Perform K-fold cross validation##
  folds.rows <- cut(1:p, breaks=nFolds, labels=FALSE)
  Y.list <- vector(mode = "list", length = nFolds)    #Folds of Y as a list
  for (f in 1:nFolds) {
    Y.list[[f]] <- Y[folds.rows==f,]
  }
  rm(Y)
  
  if (is.null(n_cores)) {n_cores <- max(detectCores() - 1, 1)}
  cl <- makeCluster(n_cores)
  clusterExport(cl, c("SYY", "Lambda", "maxK", "p", "tol.rho", "max.iter.rho", "svd.method"), envir=environment())
  clusterEvalQ(cl, {library(irlba); library(CorrConf)})
  out.parallel <- parSapply(cl=cl, Y.list, XVal_K.simrho)
  #out.parallel <- parLapply(cl=cl, Y.list, XVal_K.simrho)
  stopCluster(cl)
  
  out$LOO.XV <- 1/n/p * apply(out.parallel, 1, sum)
  #out$rho <- matrix(0, nrow=nFolds, ncol=maxK+1)
  #out$LOO.XV <- 0
  #for (i in 1:nFolds){ tmp <- out.parallel[[i]]; out$LOO.XV <- out$LOO.XV + 1/n/p*tmp$Loss; out$rho[i,] <- tmp$rho }
  out$K.hat <- out$K[which.min(out$LOO.XV)]
  if (plotit) {
    plot(out$K, out$LOO.XV, xlab="K", ylab="LOO-XV", main="Leave one out cross validation"); lines(out$K, out$LOO.XV)
    points(out$K.hat, min(out$LOO.XV), pch="x", col="red")
  }
  
  return(out)
}




######Functions to parallelize######

XVal_K.simrho <- function(Y.test) {
  p.1 <- nrow(Y.test)
  p.0 <- p - p.1
  
  SYY.0 <- 1/p.0 * (p * SYY - t(Y.test) %*% Y.test)
  train.i <- Optimize.Theta.simrho(SYY = SYY.0, Lambda = Lambda, maxK = maxK, tol.rho = tol.rho, max.iter.rho = max.iter.rho, svd.method = svd.method)
  test.loo.i <- Test.LOOXV.simrho(Y.0=Y.test, Lambda=Lambda, train=train.i)
  return(test.loo.i$Loss)   #Return only the loss. To check accuracy of code, comment this out and return test.loo.i; make sure to then use parLapply above
  #return(test.loo.i)
}

###Leave one-out cross validation###

Test.LOOXV.simrho <- function(Y.0, Lambda, train) {
  K <- train$K
  C.list <- train$C
  n <- ncol(Y.0)
  p.0 <- nrow(Y.0)
  out <- list()
  out$Loss <- rep(0, length(K))
  out$rho <- train$rho
  
  for (k in K) {
    rho.k <- out$rho[k+1]   #I will instead use the training set rho's
    V.k <- 1 + rho.k * (Lambda - 1)
    V.k <- exp(-1/n * sum(log(V.k))) * V.k   #Ensure log|V.k| = 0
    if (k == 0) {
      out$Loss[k+1] <- sum(sweep(Y.0, MARGIN = 2, 1/sqrt(V.k), "*")^2)
    } else {
      C.k <- C.list[[k+1]]
      
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

