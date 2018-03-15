\name{EstimateC_complete}
\alias{EstimateC_complete}

\title{
Estimate the n x K matrix of latent factors, C.
}
\description{
The estimate for C should be included as additional covariates in all downstream analyses.
}
\usage{
EstimateC_complete(Y, K, X = NULL, Z = NULL, B = NULL, A.ine = NULL, c.ine = NULL, A.equ = NULL, Var.0 = NULL, Cperp = NULL, rho = NULL, return.all = T, EstVariances = F, simpleDelta = F, tol.rho = 0.001, max.iter.rho = 15, return.Bhat = F, svd.method = "fast")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
A p x n data matrix; p = number of units (i.e. genes), n = sample size
}
  \item{K}{
The dimension of the latent factors, i.e. the number of columns in C
}
  \item{X}{
An n x d matrix of the covariates of interest (i.e. disease status)
}
  \item{Z}{
An n x r matrix of nuisance covariates that the user does not care about, like the intercept.
}
  \item{B}{
A list of length b of n x n positive semi-definite relatedness matrices if b > 1 or a matrix if b = 1. If the user does not supply the identity I_n or if I_n is not in the span of the other relatedness matrices, I_n will be automatically added to the list of relatedness matrices. If left unspecified, it is assumed all individuals (i.e. samples) are uncorrelated.
}
  \item{A.ine}{
A #Inequality-constraints x b matrix giving the inequality constraints on the variance multipliers. If left unspecified, it defaults to assuming all components are >= 0.
}
  \item{c.ine}{
A #Inequality-constraints length vector, where A.ine tau >= c.ine. The default is c.ine = 0, and it is highly recommended the user not define anything else.
}
  \item{A.equ}{
A #equality-constraints x b matrix giving the equality constraints on the variance multipliers. If left unspecified, it defaults to assuming there are no equality constraints. If this is specified, a using MUST define \code{Var.0} and it must be in the feasible region.
}
  \item{Var.0}{
A starting point for tau. If there are equality constraints, the user MUST specifiy a starting point within the feasible region.
}
  \item{Cperp}{
An estimate of Q.Z' C.perp. If left NULL, this is estimated with EstimateCperp. It is recommended the user leave this blank.
}
  \item{rho}{
A vestige of previous versions and should be left blank.
}
  \item{return.all}{
If \code{TRUE}, it returns all C.perp's with latent dimensions 1 up to K.
}
  \item{EstVariances}{
If \code{TRUE} and B has only ONE relatedness matrix, it estimates the variance components for every unit g = 1,...,p
}
  \item{simpleDelta}{
If \code{TRUE}, then it assumes Delta = delta^2 * I_p. Otherwise, Delta = diag(delta_1^2,...,delta_p^2). We recommend the user assign this \code{TRUE} for datasets with large p and n, i.e. methylation array data. Otherwise, the program will be prohibitively slow.
}
  \item{tol.rho}{
Each iteration of sequential PCA completes when || tau_j/||tau_j||_2 - tau_{j-1}/||tau_{j-2}||_2 ||_2 <= b * tol.rho
}
  \item{max.iter.rho}{
Maximum number of iterations j for each sequential PCA iteration
}
  \item{return.Bhat}{
If \code{TRUE}, it returns a p x d matrix Bhat, an estimate for effect of X on each unit g = 1,...,p. It is recommended the user only use this for exploratory analysis, since this is NOT the generalised least square estimate for the main effect and the variance of the estimate can only be approximated. It defaults to \code{TRUE} when samples i = 1,...,n are assumed to be uncorrelated, i.e. when B is unspecified.
}
  \item{svd.method}{
"fast": uses the package \code{irlba} for svd; any other options uses the R default.
}
}
\details{

}
\value{
A list
\item{X}{The n x d matrix of covariates of interest}
\item{Z}{The n x r matrix of nuisance covariates}
\item{Cperp}{A list of length K+1 who's entries contain the part of C orthogonal to X, with Z rotated out, when the image of C is assumed to have dimension k = 0,1,...,K}
\item{rho}{A b x (K+1) matrix or (K+1) vector containing the normalized variance components for each B_i, i.e. tau/||tau||_2}
\item{Omega.GLS}{The K x d bias-corrected estimate of CP_{Z.perp}X(X'P_{Z.perp}X)^{-1}, where P_{Z.perp} is the orthogonal project matrix onto the space orthogonal to Z}
\item{Omega.GLS.naive}{The K x d uncorrected version of Omega.GLS. This is NOT an accurate estimate and should only be used as a comparison to Omega.GLS}
\item{C}{An estimate of the n x K matrix of latent factors.}
\item{Bhat}{A p x d matrix containing the estimated effect of X on the gth unit's response, for g = 1,...,p. This is only returned when the relatedness matrices/list is \code{NULL} or when \code{return.Bhat=TRUE}.}
\item{tscores}{A p x d matrix of t-statistics for the effect of X on the gth unit's response. This is only returned when the relatedness matrices/list is \code{NULL} or when \code{return.Bhat=TRUE}.}
\item{zscores}{A p x d matrix of z-statistics for the effect of X on the gth unit's response. These should be approximately N(0,1) under the null. This is only returned when the relatedness matrices/list is \code{NULL} or when \code{return.Bhat=TRUE}.}
\item{pvalues}{A p x d matrix of p-values for each main effect. This is only returned when the relatedness matrices/list is \code{NULL} or when \code{return.Bhat=TRUE}.}
}
\references{
https://github.com/chrismckennan/CorrConf
}
\author{
Chris McKennan
}
\note{

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
https://github.com/chrismckennan/CorrConf
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (Y, K, X = NULL, Z = NULL, B = NULL, A.ine = NULL, c.ine = NULL, 
    A.equ = NULL, Var.0 = NULL, Cperp = NULL, rho = NULL, return.all = T, 
    EstVariances = F, simpleDelta = F, tol.rho = 0.001, max.iter.rho = 15, 
    return.Bhat = F, svd.method = "fast") 
{
    if (is.list(B) && length(B) > 1) {
        B <- IncludeIdent(B)
        D.ker <- CreateD.ker(A.equ)
    }
    if (is.null(X)) {
        out <- EstimateCperp(Y = Y, K = K, X = X, Z = Z, B = B, 
            simpleDelta = simpleDelta, A.ine = A.ine, c.ine = c.ine, 
            A.equ = A.equ, Var.0 = Var.0, return.all = T, tol.rho = tol.rho, 
            max.iter.rho = max.iter.rho, svd.method = svd.method)
        out$X <- X
        out$Z <- Z
        if (!is.null(Z)) {
            Q.Z <- qr.Q(qr(Z), complete = T)[, (ncol(Z) + 1):nrow(Z)]
            out$C.all <- lapply(out$C, function(C) {
                if (is.null(C)) {
                  return(C)
                }
                Q.Z \%*\% C
            })
        }
        else {
            out$C.all <- out$C
        }
        out$C <- out$C.all[[K + 1]]
        out$Cperp <- NULL
        return(out)
    }
    if (is.null(B)) {
        return.Bhat <- T
    }
    out <- list()
    out$rho <- rho
    out$Sigma.e <- NULL
    out$Sigma.b <- NULL
    out$X <- X
    out$Z <- Z
    out$Cperp <- Cperp
    if (is.null(Cperp) || (!is.null(B) && is.null(rho))) {
        out.perp <- EstimateCperp(Y = Y, K = K, X = X, Z = Z, 
            B = B, simpleDelta = simpleDelta, A.ine = A.ine, 
            c.ine = c.ine, A.equ = A.equ, Var.0 = Var.0, return.all = return.all, 
            tol.rho = tol.rho, max.iter.rho = max.iter.rho, svd.method = svd.method)
        out$Cperp <- out.perp$C
        if (return.all) {
            Cperp <- out.perp$C[[K + 1]]
        }
        else {
            Cperp <- out.perp$C
        }
        if (is.null(rho)) {
            out$rho <- out.perp$rho
            if (return.all) {
                if (is.matrix(out.perp$rho)) {
                  rho <- out.perp$rho[K + 1, ]
                }
                else {
                  rho <- out.perp$rho[K + 1]
                }
            }
            else {
                rho <- out.perp$rho
            }
        }
    }
    if (is.list(B) && length(B) == 1) {
        B <- B[[1]]
    }
    if (!is.null(Z)) {
        Q.Z <- qr.Q(qr(Z), complete = T)[, (ncol(Z) + 1):nrow(Z)]
        X <- t(Q.Z) \%*\% X
        Y <- Y \%*\% Q.Z
        if (!is.null(B)) {
            if (is.list(B)) {
                B <- lapply(B, function(x) {
                  t(Q.Z) \%*\% x \%*\% Q.Z
                })
            }
            else {
                B <- t(Q.Z) \%*\% B \%*\% Q.Z
            }
        }
    }
    p <- nrow(Y)
    n <- ncol(Y)
    d <- ncol(X)
    Q.X <- qr.Q(qr(X), complete = T)[, (d + 1):n]
    if (K == 0) {
        simpleDelta <- F
    }
    if (simpleDelta && !is.null(B)) {
        Y2 <- Y \%*\% Q.X
        if (is.list(B)) {
            out.seq <- seq.PCA.multB(Y = Y2, B = lapply(B, function(x, 
                Q.X) {
                t(Q.X) \%*\% x \%*\% Q.X
            }, Q.X = Q.X), K = K, Rho.0 = rho, A = A.ine, c = c.ine, 
                D.ker = D.ker, max.iter = 1, svd.method = svd.method)
            rho <- out.seq$Rho
            Delta.0 <- out.seq$Delta
        }
        else {
            out.seq <- seq.PCA(Y = Y2, K = K, B = t(Q.X) \%*\% 
                B \%*\% Q.X, rho.0 = rho, max.iter = 1)
            rho <- out.seq$rho
            Delta.0 <- out.seq$Delta
        }
        V <- EstimateV.complete(rho, B)
        V.tilde <- t(Q.X) \%*\% V \%*\% Q.X
        V.tilde.inv <- solve(V.tilde)
        sqrt.V.tilde <- sqrt.mat(V.tilde.inv)
        Y2 <- Y2 \%*\% sqrt.V.tilde
        if (svd.method == "fast") {
            Cperp <- sqrt(n) * cbind(svd(sqrt.mat(V.tilde) \%*\% 
                cbind(irlba(A = Y2/sqrt(Delta.0), nv = K, tol = 1/sqrt(n) * 
                  1e-04)$v[, 1:K]))$u)
        }
        else {
            Cperp <- sqrt(n) * cbind(svd(sqrt.mat(V.tilde) \%*\% 
                cbind(svd(Y2/sqrt(Delta.0), nv = K)$v[, 1:K]))$u)
        }
        Cperp <- Q.X \%*\% Cperp
        if (return.all) {
            out$Cperp[[K + 1]] <- Cperp
            if (is.list(B)) {
                out$rho[K + 1, ] <- rho
            }
            else {
                out$rho[K + 1] <- rho
            }
        }
        else {
            out$Cperp <- Cperp
            out$rho <- rho
        }
    }
    if (!simpleDelta || is.null(B)) {
        V <- EstimateV.complete(rho, B)
        if (is.null(V)) {
            V <- diag(n)
        }
        V.tilde <- t(Q.X) \%*\% V \%*\% Q.X
        V.tilde.inv <- solve(V.tilde)
        sqrt.V.tilde <- sqrt.mat(V.tilde.inv)
        Y2 <- Y \%*\% (Q.X \%*\% sqrt.V.tilde)
    }
    if (K == 0) {
        if (return.Bhat) {
            out$Bhat <- Y \%*\% solve(V, X) \%*\% solve(t(X) \%*\% 
                solve(V, X))
            out$Delta.hat <- rowSums(Y2^2)/(n - d - K)
            out$tscores <- sweep(x = out$Bhat/sqrt(out$Delta.hat), 
                MARGIN = 2, STATS = sqrt(diag(solve(t(X) \%*\% 
                  solve(V, X)))), FUN = "/", check.margin = F)
            out$zscores <- qnorm(pt(out$tscores, df = n - d - 
                K))
            out$pvalues <- 2 * pt(-abs(out$tscores), df = n - 
                d - K)
        }
        out$C <- NULL
        out$Omega.GLS <- NULL
        out$Omega.GLS.naive <- NULL
        return(out)
    }
    Cperp.reduced <- sqrt.V.tilde \%*\% t(Q.X) \%*\% Cperp
    var.mat <- solve(t(Cperp.reduced) \%*\% Cperp.reduced)
    L.hat <- Y2 \%*\% (Cperp.reduced \%*\% var.mat)
    Resids2 <- Y2 - L.hat \%*\% t(Cperp.reduced)
    Delta.hat <- rowSums(Resids2^2)/(n - d - K)
    Y1 <- Y \%*\% solve(V, X) \%*\% solve(t(X) \%*\% solve(V, X))
    out$Omega.GLS <- solve(t(L.hat/Delta.hat) \%*\% L.hat - p * 
        var.mat, t(L.hat/Delta.hat) \%*\% Y1)
    out$Omega.GLS.naive <- solve(t(L.hat/Delta.hat) \%*\% L.hat, 
        t(L.hat/Delta.hat) \%*\% Y1)
    out$C <- X \%*\% t(out$Omega.GLS) + V \%*\% Q.X \%*\% solve(t(Q.X) \%*\% 
        V \%*\% Q.X, t(Q.X) \%*\% Cperp)
    if (!is.null(Z)) {
        out$C <- Q.Z \%*\% out$C
    }
    if (return.Bhat) {
        out$Bhat <- Y1 - L.hat \%*\% out$Omega.GLS
        out$tscores <- sweep(x = out$Bhat/sqrt(Delta.hat), MARGIN = 2, 
            STATS = sqrt(diag(solve(t(X) \%*\% solve(V, X))) + 
                diag(t(out$Omega.GLS) \%*\% var.mat \%*\% out$Omega.GLS)), 
            FUN = "/", check.margin = F)
        out$zscores <- qnorm(pt(out$tscores, df = n - d - K))
        out$pvalues <- 2 * pt(-abs(out$tscores), df = n - d - 
            K)
        out$Delta.hat <- Delta.hat
    }
    if (EstVariances && is.matrix(B)) {
        inv.sqrtV.X <- sqrt.mat2((1 - rho) * diag(n - d) + rho * 
            t(Q.X) \%*\% B \%*\% Q.X)$Rinv
        Cperp.inv.X <- inv.sqrtV.X \%*\% t(Q.X) \%*\% Cperp
        Y2 <- Y \%*\% (Q.X \%*\% inv.sqrtV.X)
        L.0 <- Y2 \%*\% (Cperp.inv.X \%*\% solve(t(Cperp.inv.X) \%*\% 
            Cperp.inv.X))
        Resids <- Y2 - L.0 \%*\% t(Cperp.inv.X)
        Delta <- rowSums(Resids^2)/(n - d - K)
        rm(Resids, Y2, L.0)
        out.var <- Gene.Variances_turbo(Y = Y, Cov = cbind(X, 
            Cperp), B = B, Sigma.start = (1 - rho) * Delta, Sigma.b.start = rho * 
            Delta, tol = 1e-06)
        out$Sigma.e <- out.var$Sigma.0
        out$Sigma.b <- out.var$Sigma.b
    }
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Factor analysis }% use one of  RShowDoc("KEYWORDS")
