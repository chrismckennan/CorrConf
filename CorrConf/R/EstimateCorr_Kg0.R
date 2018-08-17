require(irlba)


#######Estimate rho and Q'C for each k = 0, 1,...,K, where Q is the orthogonal basis for ker(C')#######
Optimize.Theta.multB <- function(Y, maxK, B, Cov=NULL, A=NULL, c=NULL, D.ker=NULL, Var.0=NULL, tol.rho=1e-3, max.iter.rho=10, svd.method="fast") {
  maxK <- max(0, maxK); maxK <- round(maxK)
  if (!is.null(Cov)) {
    Q <- qr.Q(qr(Cov), complete=T)[,(ncol(Cov)+1):nrow(Cov)]
    Y <- Y %*% Q    #No more covariates
    B <- lapply(B, function(x, Q){t(Q) %*% x %*% Q}, Q=Q)
  }
  p <- nrow(Y)
  n <- ncol(Y)
  b <- length(B)
  
  out <- list()
  out$K <- 0:maxK
  out$Rho <- matrix(0, nrow=maxK+1, ncol=b)
  out$C <- vector("list", maxK+1)
  
  #K = 0#
  out.K0 <- Est.Corr.multB(Y=Y, B=B, theta.0=Var.0, simple.rho=F, A=A, c=c, D.ker=D.ker)
  out$Rho[1,] <- out.K0$Rho
  if (maxK == 0) {
    return(out)
  }
  
  #K > 0#
  Rho.0 <- out$Rho[1,]
  Delta.0 <- out.K0$Delta
  for (k in 1:maxK) {
    out.k <- seq.PCA.multB(Y=Y, B=B, K=k, Rho.0=Rho.0, Delta.0=Delta.0, A=A, c=c, D.ker=D.ker)
    Rho.0 <- out.k$Rho
    Delta.0 <- out.k$Delta
    V.0 <- CreateV(B=B, Rho=Rho.0)
    out.sqrt.V <- sqrt.mat2(V.0); sqrt.V <- out.sqrt.V$R; sqrt.Vinv <- out.sqrt.V$Rinv
    
    if (svd.method == "fast") {
      out$C[[k+1]] <- sqrt(n) * cbind(qr.Q(qr( sqrt.V %*% cbind(svd.wrapper((Y %*% sqrt.Vinv) / sqrt(Delta.0), nu=0, nv=k)$v) )))
    } else {
      out$C[[k+1]] <- sqrt(n) * cbind(qr.Q(qr( sqrt.V %*% cbind(svd.wrapper((Y %*% sqrt.Vinv) / sqrt(Delta.0), nu=0, nv=k)$v) )))
    }
    out$Rho[k+1,] <- Rho.0
  }
  return(out)
}


#######Estimate rho and Delta with sequential PCA#######
#This is only for a given K and good starting point for rho
#The convergence criterion is the Cauchy-ness of rho

seq.PCA.multB <- function(Y, B, K, Rho.0, Delta.0=NULL, A=NULL, c=NULL, D.ker=NULL, svd.method="fast", max.iter=10, tol.rho=1e-3) {
  n <- ncol(Y)
  p <- nrow(Y)
  b <- length(B)
  if (is.null(Delta.0)) {Delta.0 <- rep(1,p)}
  
  Rho.mat <- matrix(0, nrow=max.iter+1, ncol=b)
  Rho.mat[1,] <- Rho.0
  for (i in 1:max.iter) {
    V.0 <- CreateV(B, Rho.0)
    out.sqrt.V <- sqrt.mat2(V.0)
    sqrt.V <- out.sqrt.V$R; sqrt.Vinv <- out.sqrt.V$Rinv
    
    if (svd.method == "fast") {
      s.0 <- svd.wrapper((Y %*% sqrt.Vinv) / sqrt(Delta.0), nu=0, nv=K)
    } else {
      s.0 <- svd.wrapper((Y %*% sqrt.Vinv) / sqrt(Delta.0), nu=0, nv=K)
    }
    C.0 <- sqrt.V %*% s.0$v[,1:K]
    Q.C <- qr.Q(qr(C.0), complete=T)[,(K+1):n]
    
    out.rho.1 <- Est.Corr.multB(Y=Y %*% Q.C, B=lapply(B, function(x, Q.C){t(Q.C) %*% x %*% Q.C}, Q.C=Q.C), theta.0=Rho.0, A=A, c=c, D.ker=D.ker)
    Rho.1 <- out.rho.1$Rho
    if (norm(Rho.0/norm(Rho.0,type="2")-Rho.1/norm(Rho.1,type="2"), type="2") < b*tol.rho && i > 1) {
      Rho.mat[i+1,] <- Rho.1
      return(list(Rho=Rho.1, Delta=out.rho.1$Delta, all.Rho=Rho.mat[1:(i+1),], out=1))
    }
    Rho.0 <- Rho.1
    Delta.0 <- out.rho.1$Delta
    Rho.mat[i+1,] <- Rho.0
  }
  return(list(Rho=Rho.0, Delta=out.rho.1$Delta, all.Rho=Rho.mat, out=0))
}

#' Compute the square-root and inverse square-root of a symmetric, positive semi-definite matrix
#'
#' @param X An \code{n} x \code{n} symmetric psd matrix
#'
#' @return A list \item{R}{An \code{n} x \code{n} matrix; R\%*\%R = X} \item{Rinv}{An \code{n} x \code{n} matrix; Rinv\%*\%Rinv = X^\{-1\}}
#' @export
sqrt.mat2 <- function(X) {  #R^2 = X
  s <- svd.wrapper(X)
  return( list(R=sweep(s$u, 2, sqrt(s$d), "*") %*% t(s$u), Rinv=sweep(s$u, 2, 1/sqrt(s$d), "*") %*% t(s$u) ) )
}

#' Compute the square-root of a symmetric, positive semi-definite matrix
#'
#' @param X An \code{n} x \code{n} symmetric psd matrix
#'
#' @return A matrix V such that V\%*\%V = X
#' @export
sqrt.mat <- function(X) {
  s <- svd.wrapper(X)
  return(sweep(s$u, 2, sqrt(s$d), "*") %*% t(s$u))
}

#' Create a covariance matrix V from a list of matrices
#'
#' @param B A list of positive semi-definite matrices
#' @param Rho A vector of variance multipliers. length(B) = length(Rho)
#'
#' @return A covariance matrix
#' @export
CreateV <- function(B, Rho) {
  n <- ncol(B[[1]])
  V <- matrix(0, nrow=n, ncol=n)
  for (j in 1:length(B)) {
    V <- V + Rho[j]*B[[j]]
  }
  return(V)
}

##Estimate Delta when Rho and C are known
Delta.multB <- function(Y, B, Rho, Cov) {
  Q <- qr.Q(qr(Cov), complete=T)[,(ncol(Cov)+1):nrow(Cov)]
  return(Delta.K0.multB(Y = Y %*% Q, B = lapply(B, function(x, Q){t(Q) %*% x %*% Q}, Q=Q), Rho = Rho))
}
