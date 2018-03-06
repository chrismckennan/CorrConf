require(parallel)
require(irlba)

ChooseK <- function(Y, Cov=NULL, maxK=20, B=NULL, nFolds=10, simpleDelta=F, A=NULL, c=NULL, tol.rho=1e-3, max.iter.rho=15, svd.method="fast", plotit=T) {
  
  ##No random effect, samples are uncorrelated##
  if (is.null(B)) {
    return( ChooseK_NoB(Y=Y, X=Cov, maxK=maxK, nFolds=nFolds, simpleDelta=simpleDelta, max.iter.svd=3, svd.method="fast", plotit=plotit) )
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
    if (norm(B[[1]] - diag(nrow(B[[1]])), type="2") > 1e-8) {B <- c( list(diag(nrow(B[[1]]))), B )}
    if (simpleDelta) {  #Simple rho
      return( ChooseK_parallel.multB.simrho( Y=Y, X=Cov, maxK=maxK, B=B, nFolds=nFolds, A.lin=A, c.lin=c, tol.rho=tol.rho, max.iter.rho=max.iter.rho, svd.method=svd.method, plotit=plotit ) )
    } else {
      return( ChooseK_parallel.multB(Y=Y, X=Cov, maxK=maxK, B=B, nFolds=nFolds, A.lin=A, c.lin=c, tol.rho=tol.rho, max.iter.rho=max.iter.rho, svd.method=svd.method, plotit=plotit) )
    }
  }
}
