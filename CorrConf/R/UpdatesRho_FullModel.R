####This script provides all functions for seq.PCA when the delta's are allowed to be different


##Compute profile likelihood for K=0 (actually 1/(np) * PL)
#The only inputs I need are Y and rho, which can be a vector
PL.Ke0 <- function(Y, eigs.K, rho) {
  n.rho <- length(rho)
  n <- ncol(Y)
  p <- nrow(Y)
  
  if (n.rho > 1) {
    V.mat <- 1 + cbind(eigs.K - 1) %*% rbind(rho)  #n x n.rho eigenvalues of V
    ll.vec <- rep(NA, n.rho)
    
    for (i in 1:n.rho) {
      Delta.i <- 1/n * rowSums( sweep(Y, MARGIN = 2, 1/V.mat[,i], FUN="*") * Y )
      ll.vec[i] <- -1/p * sum(log(Delta.i)) - 1/n * sum(log(V.mat[,i]))
    }
    return(ll.vec - 1)
  } else {
    V <- 1 + rho * (eigs.K - 1)
    Delta <- 1/n * rowSums( sweep(Y, MARGIN = 2, 1/V, FUN="*") * Y )
    return( -1/p * sum(log(Delta)) - 1/n * sum(log(V)) - 1 )
  }
}

##Compute Delta for optimum rho when K = 0
Delta.K0 <- function(Y, eigs.K, rho) {
  n <- ncol(Y)
  V <- 1 + rho * (eigs.K - 1)
  return(1/n * rowSums( sweep(Y, MARGIN = 2, 1/V, FUN="*") * Y ))
}

##Update rho with Newton updates to MAXIMIZE the likelihood
#Repeat this until convergence
update.rho.N <- function(Y, Delta, eigs.K, rho, max.iter=1e2, tol.con=1e-4, tol.rho = 1e-4, c.wolfe=1e-3, alpha.wolfe=1/2) {
  p <- nrow(Y)
  n <- ncol(Y)
  d <- colSums((Y / Delta) * Y)    #An n-vector
  
  for (i in 1:max.iter) {
    V.i <- 1 + rho * (eigs.K - 1)
    grad.rho <- -1/n * sum((eigs.K - 1) / V.i) + 1/p/n * sum( d * (eigs.K - 1) / V.i^2 )
    if (abs(grad.rho) < tol.con || ( grad.rho > 0 && rho < tol.rho )) {   #does rho satisfy KKT conditions?
      return(list(rho=rho, n.iter=i, out=1))
    }
   
    hess.rho <- 1/n * sum((eigs.K - 1)^2 / V.i^2) - 2/p/n * sum(d * (eigs.K - 1)^2 / V.i^3)
    if (abs(hess.rho) < 1e-4) {
      step.i <- 0.1 * grad.rho / abs(grad.rho)
      step.i <- max(-rho, min(step.i, 0.99 - rho))   #Force rho + step.i to be between 0 and 0.99
    } else {
      step.i <- -grad.rho / hess.rho
      step.i <- max(-rho, min(step.i, 0.99 - rho))
    }
    
    
    ##Make sure step satisfies Wolfe conditions##
    rho.test <- rho + step.i
    V.test <- 1 + rho.test * (eigs.K - 1)
    ll.i <- -1/n * sum(log(V.i)) - 1/n/p * sum(d / V.i)
    ll.test <- -1/n * sum(log(V.test)) - 1/n/p * sum(d / V.test)
    count <- 0
    while( ll.test <= (ll.i + c.wolfe*step.i*grad.rho) && count < 10 ) {
      step.i <- alpha.wolfe * step.i
      rho.test <- rho + step.i
      V.test <- 1 + rho.test * (eigs.K - 1)
      ll.test <- -1/n * sum(log(V.test)) - 1/n/p * sum(d / V.test)
      count <- count + 1
    }
    rho <- rho.test
  }
  return(list(rho=rho, n.iter=i, out=0))
}



##Optimize rho for K = 0
#All I need is a vector of rho's to consider to get a good starting point for EM, which will refine convergence
optimize.rho.K0 <- function(Y, eigs.K, rho.start, tol=1e-6, max.iter=1e2) {
  if (length(rho.start) > 1) {
    pl.rho.start <- PL.Ke0(Y = Y, eigs.K = eigs.K, rho = rho.start)
    rho.0 <- rho.start[which.max(pl.rho.start)]    #Use this as starting point in EM updates for rho
  } else {
    rho.0 <- rho.start
  }
  
  ll.0 <- -1e16
  ll.vec <- rep(NA, max.iter+1)
  all.rhos <- rep(NA, max.iter+1)
  Newton.success <- 1
  for (i in 1:max.iter) {
    Delta.0 <- Delta.K0(Y = Y, eigs.K = eigs.K, rho=rho.0)
    out.Newton <- update.rho.N(Y = Y, Delta = Delta.0, eigs.K = eigs.K, rho = rho.0)
    if (!out.Newton$out) {
      Newton.success <- 0
    }
    rho.1 <- out.Newton$rho
    ll.1 <- PL.Ke0(Y = Y, eigs.K = eigs.K, rho = rho.1)
    ll.vec[i] <- ll.1
    all.rhos[i] <- rho.1
    if (abs(ll.1 - ll.0)/abs(ll.0) < tol) {
      return(list(rho=rho.1, out=1, n.iter=i, ll=ll.vec[1:i], all.rhos=all.rhos[1:i], out.Newton=Newton.success))
    }
    ll.0 <- ll.1
    rho.0 <- rho.1
  }
  return(list(rho=rho.1, out=0, n.iter=i))
}


