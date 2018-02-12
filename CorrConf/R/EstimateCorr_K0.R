######These functions are to be used to estimate the block correlations when K = 0######

#theta.0 is a vector of length b and B is a b-list of n x n psd matrices
#This estimates the rho's for each matrix in B and the variance multipliers, Delta
Est.Corr.multB <- function(Y, B, theta.0=NULL, max.iter=100, tol=1e-6, simple.rho = F) {
  n <- ncol(Y)
  b <- length(B)
  
  if (simple.rho) {
    if (is.null(theta.0)) {
      theta.0 <- rep(min(0.2, 1/(b+1)), b)
    }
    A <- rbind(diag(b), -diag(b), -rep(1,b))
    c <- c(rep(0,b), rep(-1,b), -1)
    out.optim <- constrOptim(theta=theta.0, f=mll.multB.simrho, grad=grad.mll.multB.simrho, ui=A, ci=c, control=list(maxit=max.iter, reltol=tol), method="BFGS", SYY=Y, B=B)
    return(list(Rho=out.optim$par, Delta=NULL, out=out.optim$convergence))
  } else {
    A <- rbind(diag(b), -diag(b), -rep(1,b))
    c <- c(rep(0,b), rep(-1,b), -1)
    
    if (is.null(theta.0)) {
      theta.0 <- rep(min(0.2, 1/(b+1)), b)
    }
    
    out.optim <- constrOptim(theta=theta.0, f=mll.multB, grad=grad.mll.multB, ui=A, ci=c, control=list(maxit=max.iter, reltol=tol), method="BFGS", Y=Y, B=B)
    Rho <- out.optim$par
    
    V <- (1-sum(Rho))*diag(n)
    for (j in 1:b) {
      V <- V + Rho[j]*B[[j]]
    }
    Delta <- 1/n * rowSums((Y %*% solve(V)) * Y)
    return(list(Rho=Rho, Delta=Delta, out=out.optim$convergence))
  }
}

##-2ll
mll.multB <- function(theta, Y, B) {
  n <- ncol(Y)
  p <- nrow(Y)
  b <- length(theta)
  
  V <- (1-sum(theta))*diag(n)
  for (j in 1:b) {
    V <- V + theta[j]*B[[j]]
  }
  s.V <- svd(V)
  Y.tilde <- Y %*% s.V$u
  YVinvYt <- rowSums( sweep(Y.tilde, 2, 1/s.V$d, FUN="*") * Y.tilde )
  return( 1/n*sum(log(s.V$d)) + 1/p*sum(log(1/n*YVinvYt)) )
}

mll.multB.simrho <- function(theta, SYY, B) {
  n <- ncol(SYY)
  b <- length(theta)
  
  V <- (1-sum(theta))*diag(n)
  for (j in 1:b) {
    V <- V + theta[j]*B[[j]]
  }
  s.V <- svd(V)
  SYY.tilde <- rowSums((t(s.V$u) %*% SYY) * t(s.V$u))
  return( 1/n * sum(log(s.V$d)) + log(1/n * sum( SYY.tilde / s.V$d )) )
}

##Gradient of -2ll
grad.mll.multB <- function(theta, Y, B) {
  n <- ncol(Y)
  p <- nrow(Y)
  b <- length(theta)
  
  V <- (1-sum(theta))*diag(n)
  for (j in 1:b) {
    V <- V + theta[j]*B[[j]]
  }
  Vinv <- solve(V)
  YVinv <- Y %*% Vinv
  YVinvYt <- rowSums(YVinv * Y)
  
  out <- rep(0, b)
  for (j in 1:b) {
    M.j <- B[[j]]; diag(M.j) <- diag(M.j) - 1
    out[j] <- 1/n * sum(M.j * Vinv) - 1/p * sum( rowSums((YVinv %*% M.j) * YVinv)/YVinvYt )
  }
  return(out)
}

grad.mll.multB.simrho <- function(theta, SYY, B) {
  n <- ncol(SYY)
  b <- length(theta)
  out <- rep(0, b)
  
  V <- (1-sum(theta))*diag(n)
  for (j in 1:b) {
    V <- V + theta[j]*B[[j]]
  }
  V.inv <- solve(V)
  SYYV.inv <- SYY %*% V.inv
  
  for (j in 1:b) {
    M.j <- B[[j]]; diag(M.j) <- diag(M.j) - 1
    M.jVinv <- M.j %*% V.inv
    out[j] <- 1/n * sum(diag(M.jVinv)) - sum( SYYV.inv * t(M.jVinv) ) / sum(diag(SYYV.inv))
  }
  return(out)
}

##Estimate Delta when Rho is known
Delta.K0.multB <- function(Y, B, Rho) {
  n <- ncol(Y)
  b <- length(Rho)
  V <- (1-sum(Rho)) * diag(n)
  for (j in 1:b) {
    V <- V + Rho[j] * B[[j]]
  }
  return(1/n * rowSums((Y %*% solve(V)) * Y))
}