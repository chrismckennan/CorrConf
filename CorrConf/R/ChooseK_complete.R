require(parallel)
require(irlba)

ChooseK <- function(Y, Cov=NULL, maxK=20, B=NULL, nFolds=10, simpleDelta=T, A.ine=NULL, c.ine=NULL, A.equ=NULL, Var.0=NULL, tol.rho=1e-3, max.iter.rho=15, svd.method="fast", plotit=T) {
  
  ##No random effect, samples are uncorrelated##
  if (is.null(B)) {
    return( ChooseK_NoB(Y=Y, X=Cov, maxK=maxK, nFolds=nFolds, simpleDelta=simpleDelta, max.iter.svd=3, svd.method=svd.method, plotit=plotit) )
  }
  
  ##One B matrix in random effect##
  if (is.matrix(B) || (is.list(B) && length(B) == 1)) {
    if (is.list(B)) {B <- B[[1]]}
    if (simpleDelta) {  #Simple rho
      return( ChooseK_parallel.simrho(Y=Y, X=Cov, maxK=maxK, B=B, nFolds=nFolds, tol.rho=tol.rho, max.iter.rho=max.iter.rho, svd.method=svd.method, plotit=plotit) )
    } else {
      return( ChooseK_parallel(Y=Y, X=Cov, maxK=maxK, B=B, nFolds=nFolds, tol.rho=tol.rho, max.iter.rho=max.iter.rho, svd.method=svd.method, plotit=plotit) )
    }
  }
  
  ##Multiple B matrices##
  #Changed notation to A.lin and c.lin here because for some reason it was not working with just A and c. I need to come back to this
  if (is.list(B) && length(B) > 1) {
    B <- IncludeIdent(B)
    D.ker <- CreateD.ker(A.equ)
    if (simpleDelta) {  #Simple rho
      return( ChooseK_parallel.multB.simrho( Y=Y, X=Cov, maxK=maxK, B=B, nFolds=nFolds, A.lin=A.ine, c.lin=c.ine, D.ker=D.ker, Var.0=Var.0, tol.rho=tol.rho, max.iter.rho=max.iter.rho, svd.method=svd.method, plotit=plotit ) )
    } else {
      return( ChooseK_parallel.multB(Y=Y, X=Cov, maxK=maxK, B=B, nFolds=nFolds, A.lin=A.ine, c.lin=c.ine, D.ker=D.ker, Var.0=Var.0, tol.rho=tol.rho, max.iter.rho=max.iter.rho, svd.method=svd.method, plotit=plotit) )
    }
  }
}


CreateD.ker <- function(A.equ) {
  if (is.null(A.equ)) { return(NULL) }
  qr.X <- qr(t(rbind(A.equ)))
  return(qr.Q(qr.X, complete=T)[,(qr.X$rank+1):ncol(A.equ)])
}

IncludeIdent <- function(B) {
  n <- nrow(B[[1]])
  mat <- cbind(c(diag(n)))
  for (i in 1:length(B)) {
    mat <- cbind(mat, c(B[[i]]))
  }
  if (qr(mat)$rank < ncol(mat)) {
    return(B)
  } else {
    return( c(list(diag(n)), B) )
  }
}
