require(irlba)

#######Estimate rho and C for each k = 0, 1,...,K#######
Optimize.Theta.multB.simrho <- function(SYY, maxK, B, Cov=NULL, A=NULL, c=NULL, D.ker=NULL, Var.0=NULL, tol.rho=1e-3, max.iter.rho=10, svd.method="fast") {
  maxK <- max(0, maxK); maxK <- round(maxK)
  if (!is.null(Cov)) {
    Q <- qr.Q(qr(Cov), complete=T)[,(ncol(Cov)+1):nrow(Cov)]
    SYY <- t(Q) %*% SYY %*% Q    #No more covariates
    B <- lapply(B, function(x, Q){t(Q) %*% x %*% Q}, Q=Q)
  }
  n <- ncol(SYY)
  b <- length(B)
  
  out <- list()
  out$K <- 0:maxK
  out$Rho <- matrix(0, nrow=maxK+1, ncol=b)
  out$C <- vector("list", maxK+1)
  
  #K = 0#
  out.K0 <- Est.Corr.multB(Y=SYY, B=B, theta.0=Var.0, simple.rho=T, A=A, c=c, D.ker=D.ker)
  out$Rho[1,] <- out.K0$Rho

  if (maxK == 0) {
    return(out)
  }
  
  #K > 0#
  Rho.0 <- out$Rho[1,]
  for (k in 1:maxK) {
    out.k <- seq.PCA.multB.simrho(SYY=SYY, B=B, K=k, Rho.0=Rho.0, A=A, c=c, D.ker=D.ker, svd.method=svd.method, max.iter=max.iter.rho, tol.rho=tol.rho)
    Rho.0 <- out.k$Rho
    V.0 <- CreateV(B=B, Rho=Rho.0)
    out.sqrt.V <- sqrt.mat2(V.0); sqrt.V <- out.sqrt.V$R; sqrt.Vinv <- out.sqrt.V$Rinv
    
    if (svd.method == "fast") {
      out$C[[k+1]] <- sqrt(n) * cbind(qr.Q(qr(sqrt.V %*% cbind(svd.wrapper(sqrt.Vinv %*% SYY %*% sqrt.Vinv, nu=0, nv=k)$v))))
    } else {
      out$C[[k+1]] <- sqrt(n) * cbind(qr.Q(qr(sqrt.V %*% cbind(svd.wrapper(sqrt.Vinv %*% SYY %*% sqrt.Vinv, nu=0, nv=k)$v))))
    }
    out$Rho[k+1,] <- Rho.0
  }
  return(out)
}


#######Estimate rho and Delta with sequential PCA#######
#This is only for a given K and good starting point for rho
#The convergence criterion is the Cauchy-ness of rho

seq.PCA.multB.simrho <- function(SYY, B, K, Rho.0, A=NULL, c=NULL, D.ker=NULL, svd.method="fast", max.iter=10, tol.rho=1e-3) {
  n <- ncol(SYY)
  b <- length(B)
  
  Rho.mat <- matrix(0, nrow=max.iter+1, ncol=b)
  Rho.mat[1,] <- Rho.0
  for (i in 1:max.iter) {
    V.0 <- CreateV(B = B, Rho = Rho.0)
    out.sqrt.V <- sqrt.mat2(V.0)
    sqrt.V <- out.sqrt.V$R; sqrt.Vinv <- out.sqrt.V$Rinv
    
    if (svd.method == "fast") {
      s.0 <- svd.wrapper(sqrt.Vinv %*% SYY %*% sqrt.Vinv, nu=0, nv=K)
    } else {
      s.0 <- svd.wrapper(sqrt.Vinv %*% SYY %*% sqrt.Vinv, nu=0, nv=K)
    }
    C.0 <- sqrt.V %*% s.0$v[,1:K]
    Q.C <- qr.Q(qr(C.0), complete=T)[,(K+1):n]
    
    tmp.1 <- t(Q.C) %*% SYY %*% Q.C; B.tmp <- lapply(B, function(x, Q.C){t(Q.C) %*% x %*% Q.C}, Q.C=Q.C)
    out.rho.1 <- Est.Corr.multB(Y=tmp.1, B=B.tmp, theta.0=Rho.0*mean(diag(tmp.1))/sum(Rho.0*sapply(B.tmp,function(x){mean(diag(x))})), simple.rho=T, A=A, c=c, D.ker=D.ker)
    rm(tmp.1,B.tmp)
    Rho.1 <- out.rho.1$Rho
    if (norm(Rho.0/norm(Rho.0,type="2")-Rho.1/norm(Rho.1,type="2"), type="2") < b*tol.rho && i > 1) {
      Rho.mat[i+1,] <- Rho.1
      return(list(Rho=Rho.1, all.Rho=Rho.mat[1:(i+1),], out=1))
    }
    Rho.0 <- Rho.1
    Rho.mat[i+1,] <- Rho.0
  }
  return(list(Rho=Rho.0, all.Rho=Rho.mat, out=0))
}

#Normalize.Rho <- 
