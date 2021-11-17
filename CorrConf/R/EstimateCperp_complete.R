###Estimate the part of C in the orthogonal complement of X###
#This returns the part of C orthogonal to X after rotating out Z, if present

require(parallel)
require(irlba)

EstimateCperp <- function(Y, K, X=NULL, Z=NULL, B=NULL, simpleDelta=F, A.ine=NULL, c.ine=NULL, A.equ=NULL, Var.0=NULL, return.all=T, tol.rho=1e-3, max.iter.rho=15, svd.method="fast") {
  if (is.list(B) && length(B) > 1) {
    B <- IncludeIdent(B)
    D.ker <- CreateD.ker(A.equ)
  }
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
    if (!is.null(X)) {X <- t(Q.Z) %*% X}
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
  
  const.y <- mean(rowMeans(Y^2))
  ######One B######
  if (is.matrix(B)) {
    if (simpleDelta) {   #Simple delta
      out.1 <- Optimize.Theta.simrho.full(SYY = 1/p/const.y * t(Y) %*% Y, X = X, B = B, maxK = K, tol.rho = tol.rho, max.iter.rho = max.iter.rho, svd.method = svd.method)
    } else {   #Multiple delta
      out.1 <- Optimize.Theta.full(Y=1/sqrt(const.y)*Y, K=K, B=B, Cov=X, tol.rho=tol.rho, max.iter.rho=max.iter.rho, svd.method=svd.method)
    }
    
    if (return.all) {
      out$C <- out.1$C
      out$K <- 0:K
      out$rho <- out.1$rho*const.y
      names(out$C) <- out$K; names(out$rho) <- out$K
    } else {
      out$C <- out.1$C[[K+1]]
      out$rho <- out.1$rho[K+1]*const.y
    }
    return(out)
  }
  
  
  ######Multiple B######
  if (is.list(B)) {
    if (simpleDelta) {   #Simple delta
      out.1 <- Optimize.Theta.multB.simrho(SYY = 1/p/const.y*t(Y)%*%Y, maxK = K, B = B, Cov = X, A=A.ine, c=c.ine, D.ker=D.ker, Var.0=Var.0, tol.rho = tol.rho, max.iter.rho = max.iter.rho, svd.method = svd.method)
    } else {   #Multiple delta
      out.1 <- Optimize.Theta.multB(Y = 1/sqrt(const.y)*Y, maxK = K, B = B, Cov = X, A=A.ine, c=c.ine, D.ker=D.ker, Var.0=Var.0, tol.rho = tol.rho, max.iter.rho = max.iter.rho, svd.method = svd.method)
    }
    if (!is.null(X) && K > 0) {
      out.1$C <- lapply(out.1$C, function(x, Q.X) {if(is.null(x)) {return(NULL)}; return(Q.X %*% x)}, Q.X=Q.X)
    }
    
    if (return.all) {
      for (k in 0:K) {out.1$Rho[k+1,] <- out.1$Rho[k+1,]/sqrt(sum(out.1$Rho[k+1,]^2))}
      out$C <- out.1$C
      out$K <- 0:K
      out$rho <- out.1$Rho
      names(out$C) <- out$K; rownames(out$rho) <- out$K
    } else {
      out$C <- out.1$C[[K+1]]
      out$rho <- out.1$Rho[K+1,]/sqrt(sum(out.1$Rho[K+1,]^2))
    }
    return(out)
  }
}

###Wrapper for SVD###

svd.wrapper <- function(x, nu=min(n,p), nv=min(n,p), symmetric=TRUE) {
  x <- as.matrix(x)
  dx <- dim(x)
  n <- dx[1L]   #x is n x p
  p <- dx[2L]
  gotit <- F
  try({svdx <- svd(x, nu, nv); gotit <- T}, silent = T)
  if(gotit) {return(svdx)}
  try({svdtx <- svd(t(x), nv, nu); gotit <- T}, silent = T)
  if (gotit) {
    temp <- svdtx$u
    svdtx$u <- svdtx$v
    svdtx$v <- temp
    return(svdtx)
  }
  if (n == p && symmetric) {
    try({temp <- eigen(x, symmetric=T); gotit <- T}, silent = T)
    if (gotit) {return(list(u=temp$vectors, v=temp$vectors, d=temp$values))}
    try({temp <- eigen(t(x), symmetric=T); gotit <- T}, silent = T)
    if (gotit) {return(list(u=temp$vectors, v=temp$vectors, d=temp$values))}
  }
  stop("SVD failed")
}
