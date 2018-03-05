######These functions are to be used to estimate the block correlations when K = 0######

#theta.0 is a vector of length b and B is a b-list of n x n psd matrices
#This estimates the rho's for each matrix in B and the variance multipliers, Delta
# A %*% theta - c >= 0
Est.Corr.multB <- function(Y, B, theta.0=NULL, A=NULL, c=NULL, max.iter=100, tol=1e-6, simple.rho = F) {
  n <- ncol(Y)
  b <- length(B)
  
  if (simple.rho) {
    if (is.null(theta.0)) {
      theta.0 <- c(1-sum(rep(min(0.2, 1/b), b-1)), rep(min(0.2, 1/b), b-1))
    }
    if (is.null(A)) {A <- diag(b)}
    if (is.null(c)) {c <- rep(1e-8,b)}
    return(0)
    out.optim <- constrOptim(theta=theta.0, f=mll.multB.simrho, grad=grad.mll.multB.simrho, ui=A, ci=c, control=list(maxit=max.iter, reltol=tol), method="BFGS", SYY=Y, B=B)
    return(list(Rho=out.optim$par, Delta=NULL, out=out.optim$convergence))
  } else {
    if (is.null(A)) {A <- diag(b)}
    if (is.null(c)) {c <- rep(1e-8,b)}
    
    if (is.null(theta.0)) {
      theta.0 <- c(1-sum(rep(min(0.2, 1/b), b-1)), rep(min(0.2, 1/b), b-1))
    }
    
    out.optim <- constrOptim(theta=theta.0, f=mll.multB, grad=grad.mll.multB, ui=A, ci=c, control=list(maxit=max.iter, reltol=tol), method="BFGS", Y=Y, B=B)
    Rho <- out.optim$par
    
    V <- CreateV(B, Rho)
    Delta <- 1/n * rowSums((Y %*% solve(V)) * Y)
    return(list(Rho=Rho, Delta=Delta, out=out.optim$convergence))
  }
}

##-2ll
mll.multB <- function(theta, Y, B) {
  n <- ncol(Y)
  p <- nrow(Y)
  b <- length(theta)
  
  V <- CreateV(B, theta)
  s.V <- svd(V)
  Y.tilde <- Y %*% s.V$u
  YVinvYt <- rowSums( sweep(Y.tilde, 2, 1/s.V$d, FUN="*") * Y.tilde )
  return( 1/n*sum(log(s.V$d)) + 1/p*sum(log(1/n*YVinvYt)) )
}

mll.multB.simrho <- function(theta, SYY, B) {
  n <- ncol(SYY)
  b <- length(theta)
  
  V <- CreateV(B, theta)
  s.V <- svd(V)
  SYY.tilde <- rowSums((t(s.V$u) %*% SYY) * t(s.V$u))
  return( 1/n * sum(log(s.V$d)) + log(1/n * sum( SYY.tilde / s.V$d )) )
}

##Gradient of -2ll
grad.mll.multB <- function(theta, Y, B) {
  n <- ncol(Y)
  p <- nrow(Y)
  b <- length(theta)
  
  V <- CreateV(B, theta)
  Vinv <- solve(V)
  YVinv <- Y %*% Vinv
  YVinvYt <- rowSums(YVinv * Y)
  
  out <- rep(0, b)
  for (j in 1:b) {
    out[j] <- 1/n * sum(B[[j]] * Vinv) - 1/p * sum( rowSums((YVinv %*% B[[j]]) * YVinv)/YVinvYt )
  }
  return(out)
}

grad.mll.multB.simrho <- function(theta, SYY, B) {
  n <- ncol(SYY)
  b <- length(theta)
  out <- rep(0, b)
  
  V <- CreateV(B, theta)
  V.inv <- solve(V)
  SYYV.inv <- SYY %*% V.inv
  
  for (j in 1:b) {
    B.jVinv <- B[[j]] %*% V.inv
    out[j] <- 1/n * sum(diag(B.jVinv)) - sum( SYYV.inv * t(B.jVinv) )/sum(diag(SYYV.inv))
  }
  return(out)
}

##Estimate Delta when Rho is known
Delta.K0.multB <- function(Y, B, Rho) {
  n <- ncol(Y)
  V <- CreateV(B, Rho)
  return(1/n * rowSums((Y %*% solve(V)) * Y))
}

Delta.K0.multB.simrho <- function(SYY, B, Rho) {
  V <- CreateV(B, Rho)
  return(1/n * sum( SYY * solve(V) ))
}
