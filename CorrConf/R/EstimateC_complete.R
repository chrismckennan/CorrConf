##This estimates the complete C, which can be used to do inference on the effect of X on Y##

library(irlba)
library(parallel)

#' Estimate \code{C} for a given value of \code{K}
#'
#' Estimate \code{C} for a given value of \code{K}. The user can also estimate Beta, the effect due to the covariates of interest by specifying \code{return.Bhat=T}.
#'
#' @param Y A \code{p} x \code{n} (number genes/methylation sites x sample size) matrix of observed expression/methylation
#' @param K The dimension of the column space of C (i.e. C is a n x K matrix). This should be estimated with ChooseK.
#' @param X A \code{n} x \code{d} design matrix of covariates of interest (i.e. disease status)
#' @param Z A \code{n} x \code{r} design matrix of nuisance covariates (i.e. known technical factors and the intercept). The default is NONE.
#' @param B A list of \code{n} x \code{n} positive semidefinite matrices that determine the correlation structure between samples. For example, these can be kinship matrices, partition matrices that group samples by blocks, diagonal matrices with 0's or 1's along the diagonal grouping samples that have the same residual variance, etc. The software automatically includes the identity if it is not in the span of the user-provided matrices. The default is to treat all samples as independent with the same residual variance.
#' @param A.ine A #inequality constaints x b matrix, where b = length(B) = num. variance components. \code{A.ine \%*\% tau >= 0}, where tau = (tau_1^2,...,tau_b^2) are the variance multipliers. The default is all variance multipliers are greater than 0.
#' @param c.ine A #inequality constraints vector, where \code{A.ine \%*\% tau - c.ine >= 0}. It is highly recommended users do NOT input anything other than \code{0}. Defaults to \code{0}.
#' @param A.equ A #equality constraints x b matrix, where \code{A.equ \%*\% tau = 0}. Defaults to no equality constraints.
#' @param Var.0 A length \code{b} starting point for tau. If A.equ is specified, this MUST be specified to ensure subsequent iterations lie in the feasible region. Otherwise, the user need not specify a starting point.
#' @param Cperp Should be left as NULL. Default is NULL
#' @param rho Should be left as NULL. Default is NULL
#' @param return.all Return variance multipliers and C for all k = 0,...,K. Default is \code{TRUE}
#' @param EstVariances If \code{TRUE} and \code{B} is a single matrix, the variance multipliers for the identity and B are estimated for each unit g=1,...,p. Default if \code{FALSE}
#' @param simpleDelta It is recommended the user not supply anything but \code{TRUE} for this argument (for now). If \code{TRUE}, the model for the residuals is MN(0, delta^2 I_p, V), where V = tau_1^2 B[[1]] + ... + tau_b^2 B[[b]]. If \code{FALSE}, Delta = diag(delta_1^2,...,delta_p^2). Defaults to \code{TRUE}.
#' @param tol.rho ICaSE (i.e. sequential PCA) terminates when || tau_j/||tau_j||_2 - tau_\{j-1\}/||tau_\{j-1\}||_2 ||_2 <= tol.rho. Default is \code{1e-3}.
#' @param max.iter.rho Maximum number of ICaSE iterations. Default is \code{15}.
#' @param return.Bhat If \code{TRUE}, the effects of interest are estimated using GLS with the model Y ~ MN(BX' + LC', Delta, V), where Delta = diag(delta.1^2,...,delta.p^2) and V are re-estimated my REML using the estimated C. Default is \code{FALSE}
#' @param svd.method Vestige of previous versions. Should not be altered by the user.
#'
#' @return A list \item{X}{The \code{n} x \code{d} model matrix} \item{Z}{The \code{n} x \code{r} model matrix of nuisance covariates, if specified} \item{C}{The estimate for the latent factors to be used in downstream analyses. The subsequent design matrix to use is [X Z C].} \item{Omega}{The estimate for the correlation between X and C} \item{Cperp}{The part of C orthogonal to \code{X}. If \code{return.all=TRUE}, it returns \code{Cperp} for all k = 1,...,K} \item{rho}{A \code{K+1} x #variance-components matrix, giving the variance multipliers for each k = 0,...,K. rho and tau are used interchangeably in this software.} \item{Delta.hat}{A p-vector of estimates for Delta, when we assume Y ~ MN_\{p x n\}(Beta X' + G Z' + L C', Delta, V) where Delta = diag(delta_1^2,...,delta_p^2). Returned only if \code{return.Bhat=TRUE}} \item{Bhat}{GLS estimate for the regression coefficients of interest, Beta, assuming Y ~ MN_\{p x n\}(Beta X' + G Z' + L C', Delta, V). If \code{B} is specified, this estimate can be improved by using a more accurate model where V is assumed gene/methylation site dependent. Returned only if \code{return.Bhat=TRUE}} \item{pvalues}{P-values for Beta. Returned only if \code{return.Bhat=TRUE}.} \item{tau.Bhat}{The tau estimated by REML under the model Y ~ MN_\{p x n\}(Beta X' + G Z' + L C', Delta, V(tau)). Returned only if \code{return.Bhat=TRUE}}
#'
#' @export
EstimateC_complete <- function(Y, K, X=NULL, Z=NULL, B=NULL, A.ine=NULL, c.ine=NULL, A.equ=NULL, Var.0=NULL, Cperp=NULL, rho=NULL, return.all=T, EstVariances=F, simpleDelta=T, tol.rho=1e-3, max.iter.rho=15, return.Bhat=F, svd.method="fast") {
  if (K <= 0) {return(0)}
  if (is.list(B) && length(B) > 1) {
    B <- IncludeIdent(B)
    D.ker <- CreateD.ker(A.equ)
  } else {
    D.ker <- NULL
  }
  if (is.null(X)) {
    out <- EstimateCperp(Y=Y, K=K, X=X, Z=Z, B=B, simpleDelta=simpleDelta, A.ine=A.ine, c.ine=c.ine, A.equ=A.equ, Var.0=Var.0, return.all=T, tol.rho=tol.rho, max.iter.rho=max.iter.rho, svd.method=svd.method)
    out$X <- X
    out$Z <- Z
    if (!is.null(Z)) {
      Q.Z <- qr.Q(qr(Z), complete = T)[,(ncol(Z)+1):nrow(Z)]
      out$C.all <- lapply(out$C, function(C) { if (is.null(C)){return(C)}; Q.Z %*% C })
    } else {
      out$C.all <- out$C
    }
    out$C <- out$C.all[[K+1]]
    out$Cperp <- NULL
    return(out)
  }
  if (is.null(B)) {return.Bhat <- T}
  
  out <- list()
  out$rho <- rho
  out$X <- X; out$Z <- Z
  out$Cperp <- Cperp
  
  if (is.null(Cperp) || (!is.null(B) && is.null(rho))) {
    out.perp <- EstimateCperp(Y=Y, K=K, X=X, Z=Z, B=B, simpleDelta=simpleDelta, A.ine=A.ine, c.ine=c.ine, A.equ=A.equ, Var.0=Var.0, return.all=return.all, tol.rho=tol.rho, max.iter.rho=max.iter.rho, svd.method=svd.method)
    out$Cperp <- out.perp$C
    if (return.all) {
      if (!is.null(B)) {
        Cperp <- out.perp$C[[K+1]]
      } else {
        Cperp <- out.perp$C
      }
    } else {
      Cperp <- out.perp$C
    }
    if (is.null(rho)) {
      out$rho <- out.perp$rho
      #out$rho <- t( apply( cbind(out$rho), 1, function(x){V.tmp <- EstimateV.complete(x,B); return( exp(-1/nrow(V.tmp) * sum(log(svd(V.tmp)$d))) * x )} ) )
      if (return.all){
        if (is.matrix(out.perp$rho)) {
          rho <- out.perp$rho[K+1,]
        } else {
          rho <- out.perp$rho[K+1]
        }
      } else {
        rho <- out.perp$rho
      }
    }
  }
  
  ##Remove the effect of Z##
  if (is.list(B) && length(B) == 1) {
    B <- B[[1]]
  }
  if (!is.null(Z)) {
    Q.Z <- qr.Q(qr(Z), complete = T)[,(ncol(Z)+1):nrow(Z)]
    X <- t(Q.Z) %*% X
    Y <- Y %*% Q.Z
    if (!is.null(B)) {
      if (is.list(B)) {
        B <- lapply(B, function(x){t(Q.Z) %*% x %*% Q.Z})
      } else {
        B <- t(Q.Z) %*% B %*% Q.Z
      }
    }
  }
  
  p <- nrow(Y)
  n <- ncol(Y)
  d <- ncol(X)
  Q.X <- qr.Q(qr(X), complete = T)[,(d+1):n]
  Y2 <- Y %*% Q.X
  
  ###Estimate Omega and C###
  if (simpleDelta) {
    V.simple <- EstimateV.complete(rho = rho, B = B)
    if (is.null(V.simple)) {V.simple <- diag(n)}
    delta.simple <- EstDelta.Simple(Y = Y2, V = t(Q.X)%*%V.simple%*%Q.X, Cov = t(Q.X)%*%Cperp)
    SCCinv.simple <- solve( t(t(Q.X)%*%Cperp) %*% solve(t(Q.X)%*%V.simple%*%Q.X, t(Q.X)%*%Cperp) )
    L.simple <- Y2 %*% solve(t(Q.X)%*%V.simple%*%Q.X, t(Q.X)%*%Cperp) %*% SCCinv.simple
    Omega.WLS.simple <- solve(t(L.simple)%*%L.simple - p*delta.simple*SCCinv.simple, t(L.simple) %*% (Y %*% solve(V.simple, X) %*% solve(t(X)%*%solve(V.simple,X))))
    out$C <- X %*% t(Omega.WLS.simple) + V.simple %*% Q.X %*% solve(t(Q.X)%*%V.simple%*%Q.X, t(Q.X)%*%Cperp)
    
    ##Return beta.hat, if prompted##
    if (return.Bhat) {
      out.Bhat <- Calc.pvalues(Y=Y, B=B, X=X, C=out$C, tau=rho, A.ine=A.ine, D.ker=D.ker)
      out$Bhat <- out.Bhat$Beta.hat
      out$pvalues <- out.Bhat$p
      out$Delta.hat <- out.Bhat$Delta.hat
      out$tau.Bhat <- out.Bhat$tau
    }
    
    if (!is.null(Z)) {out$C <- Q.Z %*% out$C}
    out$Omega <- Omega.WLS.simple
    
    return(out)
  }
  
  ##Perform 1 iteration of sequential PCA if simpleDelta is FALSE##
  if (!simpleDelta && !is.null(B)) {
    
    ##Estimate Omega.WLS with delta = diag(delta.1^2,...,delta.p^2)
    if (is.list(B)) {
      out.seq <- seq.PCA.multB(Y=Y2, B=lapply(B, function(x, Q.X){t(Q.X) %*% x %*% Q.X}, Q.X=Q.X), K=K, Rho.0=rho, A=A.ine, c=c.ine, D.ker=D.ker, max.iter=1, svd.method=svd.method)
      rho <- out.seq$Rho
      Delta.0 <- out.seq$Delta
    } else {
      out.seq <- seq.PCA(Y=Y2, K=K, B=t(Q.X) %*% B %*% Q.X, rho.0=rho, max.iter=1)
      rho <- out.seq$rho
      Delta.0 <- out.seq$Delta
    }
    V <- EstimateV.complete(rho, B)
    V.tilde <- t(Q.X) %*% V %*% Q.X; V.tilde.inv <- solve(V.tilde)
    sqrt.V.tilde <- sqrt.mat(V.tilde.inv)
    Y2 <- Y2 %*% sqrt.V.tilde
    if (svd.method=="fast") {
      Cperp <- sqrt(n) * cbind(qr.Q(qr( sqrt.mat(V.tilde) %*% cbind(svd.wrapper(Y2 / sqrt(Delta.0), nv=K, nu=0)$v[,1:K]) )))
    } else {
      Cperp <- sqrt(n) * cbind(qr.Q(qr( sqrt.mat(V.tilde) %*% cbind(svd.wrapper(Y2 / sqrt(Delta.0), nv=K, nu=0)$v[,1:K]) )))
    }
    Cperp <- Q.X %*% Cperp
    if (return.all) {
      out$Cperp[[K+1]] <- Cperp
      if (is.list(B)) {
        out$rho[K+1,] <- rho
      } else {
        out$rho[K+1] <- rho
      }
    } else {
      out$Cperp <- Cperp
      out$rho <- rho
    }
  }
  
  ######Estimate Omega with one or multiple B's######
  if (!simpleDelta || is.null(B)) {
    V <- EstimateV.complete(rho, B)
    if (is.null(V)) {V <- diag(n)}
    V.tilde <- t(Q.X) %*% V %*% Q.X; V.tilde.inv <- solve(V.tilde)
    sqrt.V.tilde <- sqrt.mat(V.tilde.inv)
    Y2 <- Y %*% (Q.X %*% sqrt.V.tilde)
  }
  
  if (K == 0) {
    if (return.Bhat) {
      out$Bhat <- Y %*% solve(V, X) %*% solve(t(X) %*% solve(V, X))
      out$Delta.hat <- rowSums(Y2^2)/(n-d-K)
      out$tscores <- sweep(x = out$Bhat / sqrt(out$Delta.hat), MARGIN = 2, STATS = sqrt(diag(solve(t(X) %*% solve(V, X)))), FUN = "/", check.margin = F)
      out$zscores <- qnorm(pt(out$tscores, df=n-d-K))
      out$pvalues <- 2*pt(-abs(out$tscores), df=n-d-K)
    }
    out$C <- NULL
    out$Omega.GLS <- NULL
    out$Omega.GLS.naive <- NULL
    return(out)
  }
  
  #I assume at least one iteration of sequential PCA has been performed#
  Cperp.reduced <- sqrt.V.tilde %*% t(Q.X) %*% Cperp; var.mat <- solve(t(Cperp.reduced) %*% Cperp.reduced)
  L.hat <- Y2 %*% (Cperp.reduced %*% var.mat)
  Resids2 <- Y2 - L.hat %*% t(Cperp.reduced); Delta.hat <- rowSums(Resids2^2)/(n-d-K)
  Y1 <- Y %*% solve(V, X) %*% solve(t(X) %*% solve(V, X))
  out$Omega.GLS <- solve(t(L.hat / Delta.hat) %*% L.hat - p*var.mat, t(L.hat / Delta.hat) %*% Y1)   #K x d
  out$Omega.GLS.naive <- solve(t(L.hat / Delta.hat) %*% L.hat, t(L.hat / Delta.hat) %*% Y1)
  
  #Estimate C#
  out$C <- X %*% t(out$Omega.GLS) + V %*% Q.X %*% solve(t(Q.X) %*% V %*% Q.X, t(Q.X) %*% Cperp)
  if (!is.null(Z)) { out$C <- Q.Z %*% out$C }
  if (return.Bhat) {
    out$Bhat <- Y1 - L.hat %*% out$Omega.GLS
    out$tscores <- sweep(x = out$Bhat / sqrt(Delta.hat), MARGIN = 2, STATS = sqrt(diag(solve(t(X) %*% solve(V, X))) + diag(t(out$Omega.GLS) %*% var.mat %*% out$Omega.GLS)), FUN = "/", check.margin = F)
    out$zscores <- qnorm(pt(out$tscores, df=n-d-K))
    out$pvalues <- 2*pt(-abs(out$tscores), df=n-d-K)
    out$Delta.hat <- Delta.hat
  }
  
  ######Estimate variances with one B######
  if (EstVariances && is.matrix(B)) {
    inv.sqrtV.X <- sqrt.mat2((1-rho)*diag(n-d) + rho*t(Q.X) %*% B %*% Q.X)$Rinv
    Cperp.inv.X <- inv.sqrtV.X %*% t(Q.X) %*% Cperp
    Y2 <- Y %*% (Q.X %*% inv.sqrtV.X)
    L.0 <- Y2 %*% (Cperp.inv.X %*% solve(t(Cperp.inv.X) %*% Cperp.inv.X))
    Resids <- Y2 - L.0 %*% t(Cperp.inv.X); Delta <- rowSums(Resids^2)/(n-d-K); rm(Resids, Y2, L.0)
    out.var <- Gene.Variances_turbo(Y=Y, Cov=cbind(X,Cperp), B=B, Sigma.start=(1-rho)*Delta, Sigma.b.start=rho*Delta, tol = 1e-6)
    out$Sigma.e <- out.var$Sigma.0
    out$Sigma.b <- out.var$Sigma.b
  }
  return(out)
}

EstimateV.complete <- function(rho, B) {
  if (is.null(B)) {
    return(NULL)
  }
  if (is.matrix(B)) {
    return( (1-rho)*diag(nrow(B)) + rho*B )
  }
  V <- matrix(0,nrow(B[[1]]),ncol(B[[1]]))
  for (j in 1:length(rho)) {
    V <- V + rho[j] * B[[j]]
  }
  return(V)
}

EstDelta.Simple <- function(Y, V, Cov=NULL) {
  if (!is.null(Cov)) {
    Q.Cov <- Compute.Q(cbind(Cov))
    V <- t(Q.Cov)%*%V%*%Q.Cov
    Y <- Y%*%Q.Cov
  }
  n <- ncol(Y)
  p <- nrow(Y)
  return(sum((Y %*% sqrt.mat2(V)$Rinv)^2) / n / p)
}

Compute.Q <- function(X) {
  qr.X <- qr(X)
  return(qr.Q(qr.X, complete = T)[,(qr.X$rank+1):nrow(X)])
}

#' Estimate beta and compute p-values under the model Y ~ MN\{ beta \%*\% X' + Gamma \%*\% Z' + L \%*\% C', diag(delta_1^2,...,delta_p^2), V(tau) \}
#' 
#' @param Y A \code{p} x \code{n} data matrix (# of genes x # of samples)
#' @param B A list of \code{n} x \code{n} positive semi-definite matrices. Default is NULL
#' @param X The \code{n} x \code{d} covariates of interest
#' @param Z A \code{n} x \code{r} matrix of nuisance covariates. Default is NULL
#' @param C The estimated \code{n} x \code{K} matrix of latent factors (see EstimateC_complete). Default is NULL
#' @param tau The starting point for tau. Default is NULL and the user need not supply one
#' @param A.ine A #inequality constaints x b matrix, where b = length(B) = num. variance components. \code{A.ine \%*\% tau >= 0}, where tau = (tau_1^2,...,tau_b^2) are the variance multipliers. The default is all variance multipliers are greater than 0.
#' @param D.ker No need to specify this
#' 
#' @return A list \item{p}{The \code{p} x \code{d} matrix of p-values} \item{Beta.hat}{The estimate for beta, a \code{p} x \code{d} matrix} \item{Delta.hat}{A \code{p}-vector of residual variances} \item{tau}{The estimate for tau}
#'
#' @export
Calc.pvalues <- function(Y, B=NULL, X, Z=NULL, C=NULL, tau=NULL, A.ine=NULL, D.ker=NULL) {
  X <- cbind(X); d <- ncol(X)
  Cov <- cbind(X,Z,C)
  Q <- Compute.Q(Cov)
  if (!is.null(B)) {
    V.params <- Est.Corr.multB(Y = Y%*%Q, B = lapply(B, function(x){t(Q)%*%x%*%Q}), theta.0 = tau, A = A.ine, D.ker=D.ker)
    Delta.hat <- V.params$Delta
  
    V.invCov <- solve(CreateV(B = B, Rho = V.params$Rho), Cov)
    SXX.inv <- solve( t(Cov) %*% V.invCov )
    Beta.hat <- Y %*% V.invCov %*% SXX.inv
    return(list(p=2*pt( -abs(Beta.hat[,1:d])/sqrt(Delta.hat)/sqrt(diag(SXX.inv)[1:d]), df = nrow(Cov)-ncol(Cov) ), tau=V.params$Rho, Beta.hat=Beta.hat[,1:d], Delta.hat=Delta.hat))
  }
  Delta.hat <- rowSums((Y%*%Q)^2)/(nrow(Cov)-ncol(Cov))
  SXX.inv <- solve( t(Cov) %*% Cov )
  Beta.hat <- Y %*% Cov %*% SXX.inv
  return(list(p=2*pt( -abs(Beta.hat[,1:d])/sqrt(Delta.hat)/sqrt(diag(SXX.inv)[1:d]), df = nrow(Cov)-ncol(Cov) ), tau=NULL, Beta.hat=Beta.hat[,1:d], Delta.hat=Delta.hat))
}
