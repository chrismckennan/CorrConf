##This estimates the complete C, which can be used to do inference on the effect of X on Y##

require(irlba)
require(parallel)

EstimateC_complete <- function(Y, K, X=NULL, Z=NULL, B=NULL, Cperp=NULL, rho=NULL, return.all=T, EstVariances=F, simpleDelta=F, tol.rho=1e-3, max.iter.rho=15, return.Bhat=F, svd.method="fast") {
  if (K == 0) {simpleDelta <- F}
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
  if (is.null(B)) {return.Bhat <- T}
  
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
  
  
  ##Perform 1 iteration of sequential PCA if simpleDelta is TRUE##
  if (simpleDelta && !is.null(B)) {
    Y2 <- Y %*% Q.X
    if (is.list(B)) {
      out.seq <- seq.PCA.multB(Y=Y2, B=lapply(B, function(x, Q.X){t(Q.X) %*% x %*% Q.X}, Q.X=Q.X), K=K, Rho.0=rho, max.iter=1)
      rho <- out.seq$Rho
      Delta.0 <- out.seq$Delta
    } else {
      out.seq <- seq.PCA(Y=Y2, K=K, B=t(Q.X) %*% B %*% Q.X, rho.0=rho, max.iter=1)
      rho <- out.seq$rho
      Delta.0 <- out.seq$Delta
    }
    V <- EstimateV.complete(rho, B)
    V.tilde <- t(Q.X) %*% V %*% Q.X; V.tilde.inv <- solve(V.tilde)
    sqrt.V.tilde <- sqrt.mat(V.tilde.inv)
    Y2 <- Y2 %*% sqrt.V.tilde
    if (svd.method=="fast") {
      Cperp <- sqrt(n) * cbind(svd(sqrt.mat(V.tilde) %*% cbind(irlba(A=Y2 / sqrt(Delta.0), nv = K, tol = 1/sqrt(n) * 1e-4)$v[,1:K]))$u)
    } else {
      Cperp <- sqrt(n) * cbind(svd(sqrt.mat(V.tilde) %*% cbind(svd(Y2 / sqrt(Delta.0), nv=K)$v[,1:K]))$u)
    }
    Cperp <- Q.X %*% Cperp
    if (return.all) {
      out$Cperp[[K+1]] <- Cperp
      if (is.list(B)) {
        out$rho[K+1,] <- rho
      } else {
        out$rho[K+1] <- rho
      }
    } else {
      out$Cperp <- Cperp
      out$rho <- rho
    }
  }
  
  ######Estimate Omega with one or multiple B's######
  if (!simpleDelta || is.null(B)) {
    V <- EstimateV.complete(rho, B)
    if (is.null(V)) {V <- diag(n)}
    V.tilde <- t(Q.X) %*% V %*% Q.X; V.tilde.inv <- solve(V.tilde)
    sqrt.V.tilde <- sqrt.mat(V.tilde.inv)
    Y2 <- Y %*% (Q.X %*% sqrt.V.tilde)
  }
  
  if (K == 0) {
    if (return.Bhat) {
      out$Bhat <- Y %*% solve(V, X) %*% solve(t(X) %*% solve(V, X))
      out$Delta.hat <- rowSums(Y2^2)/(n-d-K)
      out$tscores <- sweep(x = out$Bhat / sqrt(out$Delta.hat), MARGIN = 2, STATS = sqrt(diag(solve(t(X) %*% solve(V, X)))), FUN = "/", check.margin = F)
      out$zscores <- qnorm(pt(out$tscores, df=n-d-K))
      out$pvalues <- 2*pt(-abs(out$tscores), df=n-d-K)
    }
    out$C <- NULL
    out$Omega.GLS <- NULL
    out$Omega.GLS.naive <- NULL
    return(out)
  }
  
  #I assume at least one iteration of sequential PCA has been performed#
  Cperp.reduced <- sqrt.V.tilde %*% t(Q.X) %*% Cperp; var.mat <- solve(t(Cperp.reduced) %*% Cperp.reduced)
  L.hat <- Y2 %*% (Cperp.reduced %*% var.mat)
  Resids2 <- Y2 - L.hat %*% t(Cperp.reduced); Delta.hat <- rowSums(Resids2^2)/(n-d-K)
  Y1 <- Y %*% solve(V, X) %*% solve(t(X) %*% solve(V, X))
  out$Omega.GLS <- solve(t(L.hat / Delta.hat) %*% L.hat - p*var.mat, t(L.hat / Delta.hat) %*% Y1)   #K x d
  out$Omega.GLS.naive <- solve(t(L.hat / Delta.hat) %*% L.hat, t(L.hat / Delta.hat) %*% Y1)
  
  #Estimate C#
  out$C <- X %*% t(out$Omega.GLS) + V %*% Q.X %*% solve(t(Q.X) %*% V %*% Q.X, t(Q.X) %*% Cperp)
  if (!is.null(Z)) { out$C <- Q.Z %*% out$C }
  if (return.Bhat) {
    out$Bhat <- Y1 - L.hat %*% out$Omega.GLS
    out$tscores <- sweep(x = out$Bhat / sqrt(Delta.hat), MARGIN = 2, STATS = sqrt(diag(solve(t(X) %*% solve(V, X))) + diag(t(out$Omega.GLS) %*% var.mat %*% out$Omega.GLS)), FUN = "/", check.margin = F)
    out$zscores <- qnorm(pt(out$tscores, df=n-d-K))
    out$pvalues <- 2*pt(-abs(out$tscores), df=n-d-K)
    out$Delta.hat <- Delta.hat
  }
  
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
  if (is.null(B)) {
    return(NULL)
  }
  if (is.matrix(B)) {
    return( (1-rho)*diag(nrow(B)) + rho*B )
  }
  V <- (1-sum(rho))*diag(nrow(B[[1]]))
  for (j in 1:length(rho)) {
    V <- V + rho[j] * B[[j]]
  }
  return(V)
}
