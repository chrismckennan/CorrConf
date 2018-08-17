#require(turboEM)
require(parallel)
#########Run Turbo EM#########

#LL for turbo EM
ll.em_turbo <- function(par, y, Lambda, r) {
  Sigma <- par[2] * Lambda + par[1]
  return(-sum(log( Sigma )) - sum(y^2 / Sigma))
}

#Gradient of LL, assuming covariates have been rotated out and B is already diagonal
grad.ll <- function(y, lambda, v, v.b) {
  V <- v + v.b * lambda
  dl.dv <- sum(y^2 / V^2) - sum(1/V)
  dl.dvb <- sum(y^2 * lambda / V^2) - sum(lambda / V)
  return( c(dl.dv, dl.dvb) )
}

#Update parameters
update.vs <- function(par, y, Lambda, r) {
  v <- par[1]
  v.b <- par[2]
  y1 <- y[1:r]
  n <- length(y)
  
  Denom <- Lambda * v.b + v
  v.new <- 1/n * v^2 * sum( y^2 / Denom^2 ) + v * v.b / n * sum(Lambda / Denom)
  
  Denom <- Lambda[1:r] * v.b + v.new
  v.b.new <- 1/r * v.b^2 * sum( y1^2 * Lambda[1:r] / Denom^2 ) + 1/r * v.b * v.new * sum( 1/Denom )
  return(c(v.new, v.b.new))
}

#Turbo EM for each gene

Gene.Variances_turbo <- function(Y, Cov=NULL, B, Sigma.start, Sigma.b.start, tol = 1e-6) {
  if (!is.null(Cov)) {
    if (is.vector(B)) {
      B <- diag(B)
    }
    Q <- qr.Q(qr(Cov), complete = T)[,(ncol(Cov)+1):nrow(Cov)]
    Y <- Y %*% Q
    B <- t(Q) %*% B %*% Q 
    svd.B <- svd(B); U <- svd.B$v; Lambda <- svd.B$d
    Y <- Y %*% U
  } else {
    if (is.vector(B)) {
      Lambda <- B
    } else {
      svd.B <- svd(B); U <- svd.B$v; Lambda <- svd.B$d
      Y <- Y %*% U
    }
  }

  n <- ncol(Y)
  p <- nrow(Y)
  r <- sum(Lambda > 1e-6)   #Numerical rank of B
  
  
  Y.list <- lapply(seq_len(nrow(Y)), function(i) c(Y[i,], Sigma.start[i], Sigma.b.start[i]))
  rm(Y)
  
  cl <- makeCluster(max(detectCores() - 1, 1))
  clusterExport(cl, c("Lambda", "r", "tol"), envir=environment())
  #clusterEvalQ(cl, {source("../R/CompleteCode/OneB/IndividualVariances_parallel.R")})
  out.parallel <- parLapply(cl=cl, Y.list, Est.var)
  stopCluster(cl)
  
  out <- list()
  tmp.mat <- sapply(1:p, function(g) { tmp <- out.parallel[[g]]; c(tmp$out, tmp$Sigmas, tmp$Grad) })
  out$Out <- tmp.mat[1,]
  out$Sigma.0 <- tmp.mat[2,]; out$Sigma.b <- tmp.mat[3,]
  out$Grad <- t(tmp.mat[4:5,])
  
  return(out)
}


##Parallelized functions##
#I need Lambda, r and tol
Est.var <- function(x) {
  y <- x[1:(length(x) - 2)]
  s2.0 <- x[length(x)-1]; s2.b <- x[length(x)]
  out <- list()
  
  out.em <- turboem(par = c(s2.0, s2.b), fixptfn = update.vs, objfn=ll.em_turbo, method="squarem", Lambda=Lambda, y=y, r=r, control.run=list(tol=tol, convtype="objfn"))
  out$out <- as.numeric(!out.em$fail)
  out$Sigmas <- out.em$pars
  out$Grad <- grad.ll(y=y, lambda=Lambda, v=out$Sigmas[1], v.b=out$Sigmas[2])/length(y)
  
  return(out)
}
