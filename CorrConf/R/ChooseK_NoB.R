library(irlba)
library(parallel)

ChooseK_NoB <- function(Y, X=NULL, maxK, nFolds=10, simpleDelta=F, max.iter.svd=3, svd.method="fast", plotit=T, n_cores=NULL) {
  if (maxK < 1) {
    return(0)
  }
  
  out <- list()
  out$K <- 0:maxK
  
  if (!is.null(X)) {
    Q <- qr.Q(qr(X), complete=T)[,(ncol(X)+1):nrow(X)]
    Y <- Y %*% Q
  }
  
  p <- nrow(Y)
  n <- ncol(Y)
  nFolds <- min(nFolds, n, p)
  
  ##Permute rows and columns of Y##
  perm.rows <- order(runif(p))
  Y <- Y[perm.rows,]
  
  ##Perform K-fold cross validation##
  folds.rows <- cut(1:p, breaks=nFolds, labels=FALSE)
  
  if (is.null(n_cores)) {n_cores <- max(detectCores() - 1, 1)}
  cl <- makeCluster(n_cores)
  clusterEvalQ(cl=cl, {library(irlba); library(CorrConf)})
  if (!simpleDelta) {
    clusterExport(cl, c("Y", "maxK", "max.iter.svd", "svd.method", "folds.rows"), envir=environment())
    out.parallel <- parSapply(cl=cl, 1:nFolds, XVal_K.noB)
  } else {
    pSYY <- t(Y) %*% Y
    clusterExport(cl, c("pSYY", "maxK", "svd.method"), envir=environment())
    Y.list <- vector(mode = "list", length = nFolds)    #Folds of Y as a list
    for (f in 1:nFolds) {
        Y.list[[f]] <- Y[folds.rows==f,]
    }
    rm(Y)
    out.parallel <- parSapply(cl=cl, Y.list, XVal_K.noB.simrho)
  }
  stopCluster(cl)
  
  out <- list()
  out$K <- seq(0, maxK)
  out$LOO.XV <- rowSums(out.parallel)/(n*p)
  if (plotit) {
    plot(out$K, out$LOO.XV, xlab="K.hat", ylab="LOO X-val")
    points(out$K[which.min(out$LOO.XV)], min(out$LOO.XV), pch="x", col="red")
  }
  out$K.hat <- out$K[which.min(out$LOO.XV)]
  return(out)
}


######The function to run in parallel######

XVal_K.noB <- function(i) {

  ##Run code##
  Y.1 <- Y[folds.rows != i,];    #Train with this data (estimate C)
  Y.0 <- Y[folds.rows == i,];    #Test with these data

  train.i <- Train.noB(Y=Y.1, maxK=maxK, max.iter.svd=max.iter.svd, svd.method=svd.method)
  test.i <- Test.noB(Y.0=Y.0, train=train.i)
  return(test.i$Loss)
}

XVal_K.noB.simrho <- function(Y.i) {
  n.i <- ncol(Y.i)
  pSYY.i <- pSYY - t(Y.i) %*% Y.i
  if (svd.method == "fast") {
    s.i <- irlba(A=pSYY.i, nv=maxK, tol=1/sqrt(n.i)*1e-4)
  } else {
    s.i <- svd(pSYY.i, nv=maxK)
  }
  train.i <- list()
  train.i$C <- vector(mode="list", length=maxK)
  for (k in 1:maxK) {
    train.i$C[[k]] <- sqrt(n.i) * cbind(s.i$v[,1:k])
  }
  test.i <- Test.noB(Y.0=Y.i, train=train.i)
  return(test.i$Loss)
}

#########Test##########
#Leave one out CV
Test.noB <- function(Y.0, train) {
  C.list <- train$C
  K <- seq(0, length(C.list))
  
  out <- list()
  out$K <- K
  out$Loss <- rep(0, length(K))
  
  for (k in K) {
    if (k == 0) {
      out$Loss[k+1] <- sum(Y.0^2)
    } else {
      C.k <- C.list[[k]]
      SCC.inv <- solve(t(C.k)%*%C.k)
      L.k <- Y.0 %*% C.k %*% SCC.inv
      H.k <- apply(C.k, 1, function(x, A) {sum(x * (A %*% x))}, A=SCC.inv)
      R.k <- Y.0 - L.k %*% t(C.k)
      out$Loss[k+1] <- sum( sweep(x = R.k, MARGIN = 2, 1/(1-H.k), FUN = "*")^2 )
    }
  }
  return(out)
}


#########Train#########
Train.noB <- function(Y, maxK, max.iter.svd=3, svd.method="fast") {
  n <- ncol(Y)
  
  out <- list()
  out$C <- vector(mode="list", length=maxK)
  for (k in 1:maxK) {
    for (i in 1:max.iter.svd) {
      if (i == 1) {
        if (svd.method=="fast") {
          s.ki <- irlba(A=Y, nv=k, tol=1/sqrt(n)*1e-4)
        } else {
          s.ki <- svd(Y, nv=k)
        }
      } else {
        if (svd.method=="fast") {
          s.ki <- irlba(A=Y/sqrt(Sigma), nv=k, tol=1/sqrt(n)*1e-4)
        } else {
          s.ki <- svd(Y/sqrt(Sigma), nv=k)
        }        
      }
      C.ki <- cbind(s.ki$v[,1:k])
      
      if (i < max.iter.svd) {
        R.ki <- Y - Y %*% C.ki %*% solve(t(C.ki)%*%C.ki, t(C.ki))
        Sigma <- 1/(n-k) * rowSums(R.ki^2)
      } else {
        out$C[[k]] <- C.ki
      }
    }
  }
  return(out)
}


