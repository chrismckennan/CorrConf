require(irlba)

##This performs sequential PCA, assuming all of the delta's are the same
#This is useful when the number of variables is prohibitively large (i.e. methylation data)

#SYY is 1/p * Y'Y, and NOTHING has been rotated out
Optimize.Theta.simrho.full <- function(SYY, X=NULL, B, maxK, tol.rho=1e-3, max.iter.rho=15, svd.method="fast") {
  if (!is.null(X)) {
    Q <- qr.Q(qr(X), complete=T)[,(ncol(X)+1):nrow(X)]
    B <- t(Q) %*% B %*% Q
    SYY <- t(Q) %*% SYY %*% Q
  }
  s.B <- svd.wrapper(B)
  out <- Optimize.Theta.simrho(SYY = t(s.B$u) %*% SYY %*% s.B$u, Lambda=s.B$d, maxK=maxK, tol.rho=tol.rho, max.iter.rho=max.iter.rho, svd.method=svd.method)
  if (!is.null(X)) {
    tmp <- Q %*% s.B$u
  } else {
    tmp <- s.B$u
  }
  if (maxK > 0) {
    out$C <- lapply(out$C, function(x, A) {if(is.null(x)) {return(NULL)}; return(A %*% x)}, A=tmp)
  }
  return( out )
}

#This assumes Y has been rotated of all covariates and by the eigenvalues of Q'BQ. SYY = 1/p * U'Q' Y'Y QU
Optimize.Theta.simrho <- function(SYY, Lambda, maxK, tol.rho=1e-3, max.iter.rho=15, svd.method="fast") {
  n <- nrow(SYY)
  out <- list()
  out$K <- 0:maxK
  out$rho <- rep(0, maxK+1)
  out$C <- vector(mode = "list", length = maxK)
  
  ##K = 0##
  rho.0 <- PL.simdelta(Y=diag(SYY), Lambda=Lambda, rho=seq(0.05, 0.95, by=0.05))$rho
  out$rho[1] <- rho.0
  if (maxK == 0) {
    return(out)
  }
  
  ##Estimate when K > 0##
  for (k in 1:maxK) {
    out.k <- seq.PCA.simrho(SYY=SYY, K=k, Lambda=Lambda, rho.0=out$rho[k], tol=tol.rho, max.iter=max.iter.rho, svd.method=svd.method)
    out$rho[k+1] <- out.k$rho
    V.k <- 1 + out$rho[k+1]*(Lambda-1)
    SYY.stand <-  sweep( SYY / sqrt(V.k), MARGIN = 2, 1/sqrt(V.k), FUN = "*" )
    if (svd.method == "fast") {
      out$C[[k+1]] <- sqrt(n) * cbind(qr.Q(qr(svd.wrapper(SYY.stand, nu = 0, nv = k)$v * sqrt(V.k))))
    } else {
      out$C[[k+1]] <- sqrt(n) * cbind(qr.Q(qr(svd.wrapper(SYY.stand, nu = 0, nv = k)$v * sqrt(V.k))))
    }
  }
  return(out)
}


seq.PCA.simrho <- function(SYY, K, Lambda, rho.0, tol=1e-3, max.iter=15, svd.method="fast") {
  n <- nrow(SYY)
  rho.vec <- rep(NA, max.iter+1)
  rho.vec[1] <- rho.0
  for (i in 1:max.iter) {
    V.0 <- 1 + rho.0 * (Lambda - 1)
    SYY.tilde <- sweep(SYY / sqrt(V.0), MARGIN = 2, 1/sqrt(V.0), FUN = "*")
    if (svd.method == "fast") {
      s <- svd.wrapper(SYY.tilde, nu=0, nv=K)
    } else {
      s <- svd.wrapper(SYY.tilde, nu=0, nv=K)
    }
    Chat <- sqrt(n) * (cbind(s$v[,1:K]) * sqrt(V.0))
    Q.Chat <- qr.Q(qr(Chat), complete=T)[,(K+1):n]
    svd.QCtBQC <- svd.wrapper(t(Q.Chat * Lambda) %*% Q.Chat)
    
    #update rho#
    tmp <- Q.Chat %*% svd.QCtBQC$u
    SYY.update <- t(tmp) %*% SYY %*% tmp
    rho.1 <- Optimize.rho.simdelta(Y=diag(SYY.update), Lambda = svd.QCtBQC$d, rho.0 = rho.0)$rho
    rho.vec[i+1] <- rho.1
    if (abs(rho.1 - rho.0) < tol && i > 1) {
      return(list(rho=rho.1, n.iter=i, out=1, all.rhos=rho.vec[1:(i+1)]))
    }
    rho.0 <- rho.1
  }
  return(list(rho=rho.1, n.iter=i, out=0, all.rhos=rho.vec))
}
