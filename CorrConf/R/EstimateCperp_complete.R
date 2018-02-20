###Estimate the part of C in the orthogonal complement of X###
#This returns the part of C orthogonal to X after rotating out Z, if present

require(parallel)
require(irlba)

EstimateCperp <- function(Y, K, X=NULL, Z=NULL, B=NULL, simpleDelta=F, return.all=T, tol.rho=1e-3, max.iter.rho=15, svd.method="fast") {
  out <- list()
  out$K <- K
  out$rho <- NULL
  out$C <- NULL
  p <- nrow(Y)
  if (is.null(B) && K == 0){return(out)}
  
  if (is.list(B) && length(B) == 1) {
    B <- B[[1]]
  }
  if (!is.null(Z)) {
    Q.Z <- qr.Q(qr(Z), complete = T)[,(ncol(Z)+1):nrow(Z)]
    X <- t(Q.Z) %*% X
    Y <- Y %*% Q.Z
    if (!is.null(B)) {
      if (is.list(B)) {
        B <- lapply(B, function(x){t(Q.Z) %*% x %*% Q.Z})
      } else {
        B <- t(Q.Z) %*% B %*% Q.Z
      }
    }
  }

  if (!is.null(X)) {
    Q.X <- qr.Q(qr(X), complete=T)[,(ncol(X)+1):nrow(X)]
  }  
  
  ######No B######
  if (is.null(B)) {
    if (simpleDelta) {
      max.iter.svd <- 1
    } else {
      max.iter.svd <- 3
    }
    
    if (!is.null(X)) {
      Y <- Y %*% Q.X
    }
    n <- ncol(Y)
    
    for (i in 1:max.iter.svd) {
      if (i == 1) {
        if (svd.method=="fast") {
          s.i <- irlba(A=Y, nv=K, tol=1/sqrt(n)*1e-4)
        } else {
          s.i <- svd(Y)
        }
      } else {
        if (svd.method=="fast") {
          s.i <- irlba(A=Y/sqrt(Sigma), nv=K, tol=1/sqrt(n)*1e-4)
        } else {
          s.i <- svd(Y/sqrt(Sigma))
        }        
      }
      C.i <- cbind(s.i$v[,1:K])
      
      if (i < max.iter.svd) {
        R.i <- Y - Y %*% C.i %*% solve(t(C.i)%*%C.i, t(C.i))
        Sigma <- 1/(n-K) * rowSums(R.i^2)
      } else {
        out$C <- C.i
        if (!is.null(X)) {
          out$C <- Q.X %*% out$C
        }
      }
    }
    return(out)
  }
  
  ######One B######
  if (is.matrix(B)) {
    if (simpleDelta) {   #Simple delta
      out.1 <- Optimize.Theta.simrho.full(SYY = 1/p * t(Y) %*% Y, X = X, B = B, maxK = K, tol.rho = tol.rho, max.iter.rho = max.iter.rho, svd.method = svd.method)
    } else {   #Multiple delta
      out.1 <- Optimize.Theta.full(Y=Y, K=K, B=B, Cov=X, tol.rho=tol.rho, max.iter.rho=max.iter.rho, svd.method=svd.method)
    }
    
    if (return.all) {
      out$C <- out.1$C
      out$K <- 0:K
      out$rho <- out.1$rho
      names(out$C) <- out$K; names(out$rho) <- out$K
    } else {
      out$C <- out.1$C[[K+1]]
      out$rho <- out.1$rho[K+1]
    }
    return(out)
  }
  
  
  ######Multiple B######
  if (is.list(B)) {
    if (simpleDelta) {   #Simple delta
      out.1 <- Optimize.Theta.multB.simrho(SYY = 1/p*t(Y)%*%Y, maxK = K, B = B, Cov = X, tol.rho = tol.rho, max.iter.rho = max.iter.rho, svd.method = svd.method)
    } else {   #Multiple delta
      out.1 <- Optimize.Theta.multB(Y = Y, maxK = K, B = B, Cov = X, tol.rho = tol.rho, max.iter.rho = max.iter.rho, svd.method = svd.method)
    }
    if (!is.null(X)) {
      out.1$C <- lapply(out.1$C, function(x, A) {if(is.null(x)) {return(NULL)}; return(A %*% x)}, A=Q.X)
    }
    
    if (return.all) {
      out$C <- out.1$C
      out$K <- 0:K
      out$rho <- out.1$Rho
      names(out$C) <- out$K; rownames(out$rho) <- out$K
    } else {
      out$C <- out.1$C[[K+1]]
      out$rho <- out.1$Rho[K+1,]
    }
    return(out)
  }
}
