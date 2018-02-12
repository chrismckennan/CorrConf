##This estimates the complete C, which can be used to do inference on the effect of X on Y##

require(irlba)
require(parallel)

EstimateC_complete <- function(Y, K, X=NULL, Z=NULL, B=NULL, Cperp=NULL, rho=NULL, return.all=T, EstVariances=F, simpleDelta=F, tol.rho=1e-3, max.iter.rho=15, svd.method="fast") {
  if (is.null(X)) {
    out <- EstimateCperp(Y=Y, K=K, X=X, Z=Z, B=B, simpleDelta=simpleDelta, return.all=T, tol.rho=tol.rho, max.iter.rho=max.iter.rho, svd.method=svd.method)
    out$X <- X
    out$Z <- Z
    if (!is.null(Z)) {
      Q.Z <- qr.Q(qr(Z), complete = T)[,(ncol(Z)+1):nrow(Z)]
      out$C.all <- lapply(out$C, function(C) { if (is.null(C)){return(C)}; Q.Z %*% C })
    } else {
      out$C.all <- out$C
    }
    out$C <- out$C.all[[K+1]]
    out$Cperp <- NULL
    return(out)
  }
  
  out <- list()
  out$rho <- rho
  out$Sigma.e <- NULL; out$Sigma.b <- NULL
  out$X <- X; out$Z <- Z
  out$Cperp <- Cperp
  
  if (is.null(Cperp) || (!is.null(B) && is.null(rho))) {
    out.perp <- EstimateCperp(Y=Y, K=K, X=X, Z=Z, B=B, simpleDelta=simpleDelta, return.all=return.all, tol.rho=tol.rho, max.iter.rho=max.iter.rho, svd.method=svd.method)
    out$Cperp <- out.perp$C
    if (return.all) {
      Cperp <- out.perp$C[[K+1]]
    } else {
      Cperp <- out.perp$C
    }
    if (is.null(rho)) {
      out$rho <- out.perp$rho
      if (return.all){
        if (is.matrix(out.perp$rho)) {
          rho <- out.perp$rho[K+1,]
        } else {
          rho <- out.perp$rho[K+1]
        }
      } else {
        rho <- out.perp$rho
      }
    }
  }
  
  ##Remove the effect of Z##
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
  
  p <- nrow(Y)
  n <- ncol(Y)
  d <- ncol(X)
  Q.X <- qr.Q(qr(X), complete = T)[,(d+1):n]
  
  ######Estimate Omega with one or multiple B's######
  V <- EstimateV.complete(rho, B)
  V.tilde <- t(Q.X) %*% V %*% Q.X; V.tilde.inv <- solve(V.tilde)
  
  #I assume at least one iteration of sequential PCA has been performed#
  sqrt.V.tilde <- sqrt.mat(V.tilde.inv)
  Y2 <- Y %*% Q.X %*% sqrt.V.tilde
  Cperp.reduced <- sqrt.V.tilde %*% t(Q.X) %*% Cperp; var.mat <- solve(t(Cperp.reduced) %*% Cperp.reduced)
  L.hat <- Y2 %*% Cperp.reduced %*% var.mat
  Resids2 <- Y2 - L.hat %*% t(Cperp.reduced); Delta.hat <- rowSums(Resids2^2)/(n-d-K)
  Y1 <- Y %*% solve(V, X) %*% solve(t(X) %*% solve(V, X))
  out$Omega.GLS <- solve(t(L.hat / Delta.hat) %*% L.hat - p*var.mat, t(L.hat / Delta.hat) %*% Y1)   #K x d
  
  #Estimate C#
  out$C <- X %*% t(out$Omega.GLS) + V %*% Q.X %*% solve(t(Q.X) %*% V %*% Q.X, t(Q.X) %*% Cperp)
  if (!is.null(Z)) { out$C <- Q.Z %*% out$C }
  
  ######Estimate variances with one B######
  if (EstVariances && is.matrix(B)) {
    inv.sqrtV.X <- sqrt.mat2((1-rho)*diag(n-d) + rho*t(Q.X) %*% B %*% Q.X)$Rinv
    Cperp.inv.X <- inv.sqrtV.X %*% t(Q.X) %*% Cperp
    Y2 <- Y %*% (Q.X %*% inv.sqrtV.X)
    L.0 <- Y2 %*% (Cperp.inv.X %*% solve(t(Cperp.inv.X) %*% Cperp.inv.X))
    Resids <- Y2 - L.0 %*% t(Cperp.inv.X); Delta <- rowSums(Resids^2)/(n-d-K); rm(Resids, Y2, L.0)
    out.var <- Gene.Variances_turbo(Y=Y, Cov=cbind(X,Cperp), B=B, Sigma.start=(1-rho)*Delta, Sigma.b.start=rho*Delta, tol = 1e-6)
    out$Sigma.e <- out.var$Sigma.0
    out$Sigma.b <- out.var$Sigma.b
  }
  return(out)
}

EstimateV.complete <- function(rho, B) {
  if (is.matrix(B)) {
    return( (1-rho)*diag(nrow(B)) + rho*B )
  }
  V <- (1-sum(rho))*diag(nrow(B[[1]]))
  for (j in 1:length(rho)) {
    V <- V + rho[j] * B[[j]]
  }
  return(V)
}