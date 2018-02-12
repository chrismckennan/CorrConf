######This optimizes rho and Delta when Delta = delta * I_p######


####K = 0####
##Get estimate for delta when rho is known
est.simdelta <- function(Y, Lambda, rho) {
  if (is.vector(Y)) {
    S <- Y
    n <- length(S)
  } else {
    n <- ncol(Y)
    p <- nrow(Y)
    S <- 1/p * colSums(Y * Y)
  }
  return( 1/n * sum(S / (1 + rho * (Lambda - 1))) )
}

##Profile likelihood for rho
PL.simdelta <- function(Y, Lambda, rho) {
  if (is.vector(Y)) {
    S <- Y
    n <- length(S)
  } else {
    n <- ncol(Y)
    p <- nrow(Y)
    S <- 1/p * colSums(Y * Y)
  }
  if (length(rho) == 1) {
    V <- 1 + rho * (Lambda - 1)
    return(-n * log(sum(S / V)) - sum(log(V)))
  } else {
    ll.vec <- rep(NA, length(rho))
    
    for (i in 1:length(rho)) {
      V.i <- 1 + rho[i] * (Lambda - 1)
      ll.vec[i] <- -n * log(sum(S / V.i)) - sum(log(V.i))
    }
    rho.hat <- rho[which.max(ll.vec)]
    return(list(ll=ll.vec, rho=rho.hat, delta=est.simdelta(Y=S, Lambda=Lambda, rho=rho.hat)))
  }
}

##Newton algorithm to find rho.hat with profile likelihood
Optimize.rho.simdelta <- function(Y, Lambda, rho.0, tol=1e-5, max.iter=100, tol.rho=1e-2, alpha=0.5) {
  if (length(rho.0) > 1) {
    rho.0 <- PL.simdelta(Y=Y, Lambda=Lambda, rho=rho.0)$rho
  }
  if (is.vector(Y)) {
    S <- Y
    n <- length(S)
  } else {
    n <- ncol(Y)
    p <- nrow(Y)
    S <- 1/p * colSums(Y * Y)
  }
  for (i in 1:max.iter) {
    V.0 <- 1 + rho.0 * (Lambda - 1)
    delta.0 <- 1/n * sum(S / V.0)
    grad.0 <- 1/delta.0 * sum( S*(Lambda-1) / V.0^2 ) - sum((Lambda-1) / V.0)
    if (abs(grad.0) < n * tol || (rho.0 < tol.rho && grad.0 < 0) || (rho.0 >= 0.99 && grad.0 > 0)) {
      return(list(rho=rho.0, delta2=est.simdelta(Y=Y, Lambda=Lambda, rho=rho.0), out=1, iter=i))
    }
    hess.0 <- -2/delta.0 * sum(S * (Lambda-1)^2 / V.0^3) + 1/delta.0^2/n * sum( S*(Lambda-1) / V.0^2 )^2 + sum((Lambda-1)^2 / V.0^2)
    dir <- -grad.0/hess.0
    while (rho.0 + dir > 0.999 || rho.0 + dir < tol.rho/10) {
      dir <- alpha * dir
    }
    rho.0 <- rho.0 + dir
  }
  return(list(rho=rho.0, delta2=est.simdelta(Y=S, Lambda=Lambda, rho=rho.0), out=0, iter=i))
}