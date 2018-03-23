###This script contains functions to optimize rho using sequential PCA
#It uses the package irlba to compute the partial SVD quickly
require(irlba)

#####Estimate the projection of C onto the orthogonal complement of Cov

Optimize.Theta.full <- function(Y, K, B, Cov=NULL, tol.rho=1e-3, max.iter.rho=15, svd.method="fast") {
  if (!is.null(Cov)) {
    Q.Cov <- qr.Q(qr(Cov), complete=T)[,(ncol(Cov)+1):nrow(Cov)]
    Y <- Y %*% Q.Cov
    B <- t(Q.Cov) %*% B %*% Q.Cov
  }
  s.B <- svd.wrapper(B)
  out <- Optimize.Theta(Y=Y %*% s.B$u, K=K, B=s.B$d, tol.rho=tol.rho, max.iter.rho=max.iter.rho, svd.method=svd.method)
  if (!is.null(Cov)) {
    tmp <- Q.Cov %*% s.B$u
  } else {
    tmp <- s.B$u
  }
  out$C <- lapply(out$C, function(x, A) {if(is.null(x)) {return(NULL)}; return(A %*% x)}, A=tmp)
  return(out)
}


#####Find rho and C for all latent dimension = 0:K.
#It returns Q'C for each value of K if prompted, where Q is an orthonormal basis for ker(Cov'). Otherwise, it only returns rho

Optimize.Theta <- function(Y, K, B, Cov=NULL, tol.rho=1e-3, max.iter.rho=15, svd.method="fast", return.C=T) {
  if (!is.null(Cov)) {
    if (is.vector(B)) {
      B <- diag(B, nrow=length(B), ncol=length(B))
    }
    Q <- qr.Q(qr(Cov), complete=T)[,(ncol(Cov)+1):nrow(Cov)]
    Y <- Y %*% Q    #No more covariates
    B <- t(Q) %*% B %*% Q
    svd.B <- svd.wrapper(B)
    Lambda <- as.vector(svd.B$d)
    U <- svd.B$u
    Y <- Y %*% U      #Y now has uncorrelated columns
  } else {
    if (is.vector(B)) {
      Lambda <- B
    } else {
      svd.B <- svd.wrapper(B)
      Lambda <- as.vector(svd.B$d)
      U <- svd.B$u
      Y <- Y %*% U    #Y now has uncorrelated columns   
    }
  }
  n <- ncol(Y)
  
  out <- list()
  out$rho <- rep(0, length=K+1)
  out$K <- 0:K
  if (return.C) {
    out$C <- list(NULL)   #C = NULL when K = 0
  }
  
  ##Estimate rho when K = 0##
  rho.vec <- seq(0, 0.95, by=0.05)
  out.PL.K0 <- PL.Ke0(Y=Y, eigs.K=Lambda, rho = rho.vec)
  rho.0 <- rho.vec[which.max(out.PL.K0)]
  Delta.0 <- Delta.K0(Y=Y, eigs.K=Lambda, rho=rho.0)
  out$rho[1] <- rho.0
  if (K == 0) {
    return(out)
  }
  
  ##Estimate rho when K > 0##
  for (k in 1:K) {
    out.k <- seq.PCA(Y=Y, K=k, B=Lambda, Cov=NULL, Delta.0=Delta.0, rho.0=out$rho[k], tol=tol.rho, max.iter=max.iter.rho, svd.method=svd.method)
    out$rho[k+1] <- out.k$rho
    Delta.0 <- out.k$Delta
    Y.stand <- sweep( Y / sqrt(Delta.0), MARGIN = 2, 1/sqrt(1 + out$rho[k+1]*(Lambda-1)), FUN = "*" )
    if (return.C) {
      if (is.vector(B)) {
        if (svd.method == "fast") {
          out$C[[k+1]] <- sqrt(n) * cbind(qr.Q(qr( svd.wrapper(Y.stand, nu=0, nv=k)$v[,1:k] * sqrt(1 + out$rho[k+1]*(Lambda-1)) )))   #Q'\hat{C}
        } else {
          out$C[[k+1]] <- sqrt(n) * cbind(qr.Q(qr( svd.wrapper(Y.stand, nu=0, nv=k)$v[,1:k] * sqrt(1 + out$rho[k+1]*(Lambda-1)) )))   #Q'\hat{C}
        }
      } else {
        if (svd.method == "fast") {
          out$C[[k+1]] <- sqrt(n) * U %*% qr.Q(qr( svd.wrapper(Y.stand, nu=0, nv=k)$v[,1:k] * sqrt(1 + out$rho[k+1]*(Lambda-1)) ))   #Q'\hat{C}
        } else {
          out$C[[k+1]] <- sqrt(n) * U %*% qr.Q(qr( svd.wrapper(Y.stand, nu=0, nv=k)$v[,1:k] * sqrt(1 + out$rho[k+1]*(Lambda-1)) ))   #Q'\hat{C}
        }
      }
      
    }
  }
  return(out)
}

#####Optimize for rho and Delta for fixed K#####
#This assumes V = sigma2 * I_n + sigma2_b * B
#If B is a vector and there are no covariates, it assumes you have already rotated Y by the eigenvectors of B

seq.PCA <- function(Y, K, B, Cov=NULL, Delta.0=NULL, rho.0, tol=1e-3, max.iter=15, svd.method="fast") {
  
  if (!is.null(Cov)) {
    if (is.vector(B)) {
      B <- diag(B, nrow=length(B), ncol=length(B))
    }
    Q.cov <- qr.Q(qr(Cov), complete=T)[,(ncol(Cov)+1):nrow(Cov)]
    Y <- Y %*% Q.cov
    B <- t(Q.cov) %*% B %*% Q.cov
    svd.B <- svd.wrapper(B)
    U.cov <- svd.B$u
    Lambda.cov <- svd.B$d
    Y.bar <- Y %*% U.cov
  } else {
    if (is.vector(B)) {
      Y.bar <- Y
      Lambda.cov <- B
    } else {
      svd.B <- svd.wrapper(B)
      U.cov <- svd.B$u
      Lambda.cov <- svd.B$d
      Y.bar <- Y %*% U.cov    #I do everything with respect to Y.bar (i.e. Y.bar is the new Y)
    }
  }
  rm(Y)
  n <- ncol(Y.bar)
  p <- nrow(Y.bar)
  
  if (is.null(Delta.0)) {
    Delta.0 <- rep(1,p)
  }
  
  rho.vec <- rep(NA, max.iter+1)
  rho.vec[1] <- rho.0
  for (i in 1:max.iter) {
    V.0 <- 1 + rho.0 * (Lambda.cov - 1)
    Y.tilde.0 <- sweep(Y.bar / sqrt(Delta.0), MARGIN = 2, 1/sqrt(V.0), FUN="*")
    if (svd.method == "fast") {
      s.0 <- svd.wrapper(Y.tilde.0, nu=0, nv=K)
    } else {
      s.0 <- svd.wrapper(Y.tilde.0, nu=0, nv=K)
    }
    Chat <- sqrt(n) * (cbind(s.0$v[,1:K]) * sqrt(V.0))
    Q.Chat <- qr.Q(qr(Chat), complete=T)[,(K+1):n]
    svd.QCtBQC <- svd.wrapper(t(Q.Chat * Lambda.cov) %*% Q.Chat)
    
    #Update rho#
    Y.update <- Y.bar %*% Q.Chat %*% svd.QCtBQC$u
    rho.1 <- optimize.rho.K0(Y=Y.update, eigs.K = svd.QCtBQC$d, rho.start = rho.0)$rho
    rho.vec[i+1] <- rho.1
    Delta.0 <- Delta.K0(Y=Y.update, eigs.K = svd.QCtBQC$d, rho = rho.1)
    if (abs(rho.1 - rho.0) < tol && i > 1) {
      return(list(rho=rho.1, Delta=Delta.0, n.iter=i, out=1, all.rhos=rho.vec[1:(i+1)]))
    }
    rho.0 <- rho.1
  }
  
  return(list(rho=rho.1, Delta=Delta.0, n.iter=i, out=0, all.rhos=rho.vec))
}


