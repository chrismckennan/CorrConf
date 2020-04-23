library(parallel)
library(irlba)

#' Estimate the latent factor dimension, \code{K}, in high dimensional data
#'
#' Estimates \code{K} with a modified cross-validation that accounts for correlation between samples, if present. Of course, this method is still valid if samples are independent. The code is parallelized and runs on a user-specified number of cores.
#'
#' @param Y A \code{p} x \code{n} (number genes/methylation sites x sample size) matrix of observed expression/methylation
#' @param Cov A \code{n} x (\code{d+r}) (sample size x total number of covariates) matrix, where \code{d} = number of covariates of interest and \code{r} = number of nuisance covariates (like the intercept and/or other observed technical factors). The default is no additional covariates.
#' @param maxK The maximum number of latent factors to consider. The default is \code{20}
#' @param B An \code{n} x \code{n} positive semidefinite matrix or list of \code{n} x \code{n} positive semidefinite matrices that determine the correlation structure between samples. For example, these can be kinship matrices, partition matrices that group samples by blocks, diagonal matrices with 0's or 1's along the diagonal grouping samples that have the same residual variance, etc. The software includes the identity if it is not in the span of the user-provided matrices. The default is to treat all samples as independent with the same residual variance.
#' @param nFolds The number of gene/methylation folds to consider for cross validation. The fewer the folds, the faster the code runs and should not be set smaller than \code{2}. The smaller \code{p} is, the fewer number of folds the user should use. We recommend using no more than \code{5} for gene expression data (~1e4 genes) and no more than \code{10} for DNA methylation data (~4e5 - 8e5 CpG cites). The default is \code{5} folds.
#' @param simpleDelta The user should not supply anything but \code{TRUE} for this argument (for now). If \code{TRUE}, the model for the residuals is MN(0, delta^2 I_p, V), where V = tau_1^2 B[[1]] + ... + tau_b^2 B[[b]]. If \code{FALSE}, Delta = diag(delta_1^2,...,delta_p^2). Defaults to \code{TRUE}.
#' @param A.ine A #inequality constaints x b matrix, where b = length(B) = num. variance components. \code{A.ine \%*\% tau >= 0}, where tau = (tau_1^2,...,tau_b^2) are the variance multipliers. The default is all variance multipliers are greater than 0.
#' @param c.ine A #inequality constraints vector, where \code{A.ine \%*\% tau - c.ine >= 0}. It is highly recommended users do NOT input anything other than \code{0}. Defaults to \code{0}.
#' @param A.equ A #equality constraints x b matrix, where \code{A.equ \%*\% tau = 0}. Defaults to no equality constraints.
#' @param Var.0 A length b starting point for tau. If A.equ is specified, this MUST be specified to ensure subsequent iterations lie in the feasible region. Otherwise, the user need not specify a starting point.
#' @param tol.rho ICaSE (i.e. sequential PCA) terminates when || tau_j/||tau_j||_2 - tau_\{j-1\}/||tau_\{j-1\}||_2 ||_2 <= tol.rho. Default is \code{1e-3}.
#' @param max.iter.rho Maximum number of ICaSE iterations. Default is \code{15}.
#' @param svd.method Vestige of previous versions. Should not be altered by the user.
#' @param plotit If \code{TRUE}, plots the leave one out cross validation (LOO XV) plot. Defaults to \code{TRUE}.
#' @param n_cores The number of cores to use. The default is #available cores - 1.
#'
#' @return A list \item{K}{A vector of latent dimensions considered} \item{LOO.XV}{A vector of the average leave-one-out cross validation for each \code{K} considered} \item{K.hat}{The \code{K} that gives the minimum leave-one-out cross validation}
#' @export
ChooseK <- function(Y, Cov=NULL, maxK=20, B=NULL, nFolds=10, simpleDelta=T, A.ine=NULL, c.ine=NULL, A.equ=NULL, Var.0=NULL, tol.rho=1e-3, max.iter.rho=15, svd.method="slow", plotit=T, n_cores=NULL) {
  
  ##No random effect, samples are uncorrelated##
  if (is.null(B)) {
    return( ChooseK_NoB(Y=Y, X=Cov, maxK=maxK, nFolds=nFolds, simpleDelta=simpleDelta, max.iter.svd=3, svd.method=svd.method, plotit=plotit, n_cores=n_cores) )
  }
  
  ##One B matrix in random effect##
  if (is.matrix(B) || (is.list(B) && length(B) == 1)) {
    if (is.list(B)) {B <- B[[1]]}
    if (simpleDelta) {  #Simple rho
      return( ChooseK_parallel.simrho(Y=Y, X=Cov, maxK=maxK, B=B, nFolds=nFolds, tol.rho=tol.rho, max.iter.rho=max.iter.rho, svd.method=svd.method, plotit=plotit, n_cores=n_cores) )
    } else {
      return( ChooseK_parallel(Y=Y, X=Cov, maxK=maxK, B=B, nFolds=nFolds, tol.rho=tol.rho, max.iter.rho=max.iter.rho, svd.method=svd.method, plotit=plotit, n_cores=n_cores) )
    }
  }
  
  ##Multiple B matrices##
  #Changed notation to A.lin and c.lin here because for some reason it was not working with just A and c. I need to come back to this
  if (is.list(B) && length(B) > 1) {
    B <- IncludeIdent(B)
    D.ker <- CreateD.ker(A.equ)
    if (simpleDelta) {  #Simple rho
      return( ChooseK_parallel.multB.simrho( Y=Y, X=Cov, maxK=maxK, B=B, nFolds=nFolds, A.lin=A.ine, c.lin=c.ine, D.ker=D.ker, Var.0=Var.0, tol.rho=tol.rho, max.iter.rho=max.iter.rho, svd.method=svd.method, plotit=plotit, n_cores=n_cores ) )
    } else {
      return( ChooseK_parallel.multB(Y=Y, X=Cov, maxK=maxK, B=B, nFolds=nFolds, A.lin=A.ine, c.lin=c.ine, D.ker=D.ker, Var.0=Var.0, tol.rho=tol.rho, max.iter.rho=max.iter.rho, svd.method=svd.method, plotit=plotit, n_cores=n_cores) )
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
