\name{EstimateCperp}
\alias{EstimateCperp}

\title{
Estimate C.perp, the part of C orthogonal to X.
}
\description{
Estimate C.perp, the part of C orthogonal to X. If the n x r matrix of nuisance covariates, Z, are specified, it returns the (n-r) x K matrix Q.Z' P_{X.perp} C.hat, where P_{X.perp} is the orthogonal projection onto the space orthogonal to X and Q.Z is the n x (n-r) matrix whose columns form an orthonormal basis for ker(Z').
}
\usage{
EstimateCperp(Y, K, X = NULL, Z = NULL, B = NULL, simpleDelta = F, A.ine = NULL, c.ine = NULL, A.equ = NULL, Var.0 = NULL, return.all = T, tol.rho = 0.001, max.iter.rho = 15, svd.method = "fast")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
A p x n data matrix; p = number of units (i.e. genes), n = sample size
}
  \item{K}{
The dimension of the latent factors, i.e. the number of columns in C = number of columns in C.perp
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
  \item{simpleDelta}{
If \code{TRUE}, then it assumes Delta = delta^2 * I_p. Otherwise, Delta = diag(delta_1^2,...,delta_p^2). We recommend the user assign this \code{TRUE} for datasets with large p and n, i.e. methylation array data. Otherwise, the program will be prohibitively slow.
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
  \item{return.all}{
If \code{TRUE}, it returns all C.perp's with latent dimensions 1 up to K.
}
  \item{tol.rho}{
Each iteration of sequential PCA completes when || tau_j/||tau_j||_2 - tau_{j-1}/||tau_{j-2}||_2 ||_2 <= b * tol.rho
}
  \item{max.iter.rho}{
Maximum number of iterations j for each sequential PCA iteration
}
  \item{svd.method}{
"fast": uses the package \code{irlba} for svd; any other options uses the R default.
}
}
\details{
It is expected that most users will not use this function, as it only estimates the part of C orthogonal to X. To estimate all of C, use EstimateC_complete.
}
\value{
A list
\item{K}{The latent dimensions considered, k = 0,1,...,K}
\item{C}{A list of length K+1 who's entries contain the part of C orthogonal to X, with Z rotated out, when the image of C is assumed to have dimension k = 0,1,...,K}
\item{rho}{A b x (K+1) matrix or (K+1) vector containing the normalized variance components for each B_i}
}
\references{
See https://github.com/chrismckennan/CorrConf for more details.
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
function (Y, K, X = NULL, Z = NULL, B = NULL, simpleDelta = F, 
    A.ine = NULL, c.ine = NULL, A.equ = NULL, Var.0 = NULL, return.all = T, 
    tol.rho = 0.001, max.iter.rho = 15, svd.method = "fast") 
{
    if (is.list(B) && length(B) > 1) {
        B <- IncludeIdent(B)
        D.ker <- CreateD.ker(A.equ)
    }
    out <- list()
    out$K <- K
    out$rho <- NULL
    out$C <- NULL
    p <- nrow(Y)
    if (is.null(B) && K == 0) {
        return(out)
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
    if (!is.null(X)) {
        Q.X <- qr.Q(qr(X), complete = T)[, (ncol(X) + 1):nrow(X)]
    }
    if (is.null(B)) {
        if (simpleDelta) {
            max.iter.svd <- 1
        }
        else {
            max.iter.svd <- 3
        }
        if (!is.null(X)) {
            Y <- Y \%*\% Q.X
        }
        n <- ncol(Y)
        for (i in 1:max.iter.svd) {
            if (i == 1) {
                if (svd.method == "fast") {
                  s.i <- irlba(A = Y, nv = K, tol = 1/sqrt(n) * 
                    1e-04)
                }
                else {
                  s.i <- svd(Y)
                }
            }
            else {
                if (svd.method == "fast") {
                  s.i <- irlba(A = Y/sqrt(Sigma), nv = K, tol = 1/sqrt(n) * 
                    1e-04)
                }
                else {
                  s.i <- svd(Y/sqrt(Sigma))
                }
            }
            C.i <- cbind(s.i$v[, 1:K])
            if (i < max.iter.svd) {
                R.i <- Y - Y \%*\% C.i \%*\% solve(t(C.i) \%*\% C.i, 
                  t(C.i))
                Sigma <- 1/(n - K) * rowSums(R.i^2)
            }
            else {
                out$C <- C.i
                if (!is.null(X)) {
                  out$C <- Q.X \%*\% out$C
                }
            }
        }
        return(out)
    }
    if (is.matrix(B)) {
        if (simpleDelta) {
            out.1 <- Optimize.Theta.simrho.full(SYY = 1/p * t(Y) \%*\% 
                Y, X = X, B = B, maxK = K, tol.rho = tol.rho, 
                max.iter.rho = max.iter.rho, svd.method = svd.method)
        }
        else {
            out.1 <- Optimize.Theta.full(Y = Y, K = K, B = B, 
                Cov = X, tol.rho = tol.rho, max.iter.rho = max.iter.rho, 
                svd.method = svd.method)
        }
        if (return.all) {
            out$C <- out.1$C
            out$K <- 0:K
            out$rho <- out.1$rho
            names(out$C) <- out$K
            names(out$rho) <- out$K
        }
        else {
            out$C <- out.1$C[[K + 1]]
            out$rho <- out.1$rho[K + 1]
        }
        return(out)
    }
    if (is.list(B)) {
        if (simpleDelta) {
            out.1 <- Optimize.Theta.multB.simrho(SYY = 1/p * 
                t(Y) \%*\% Y, maxK = K, B = B, Cov = X, A = A.ine, 
                c = c.ine, D.ker = D.ker, Var.0 = Var.0, tol.rho = tol.rho, 
                max.iter.rho = max.iter.rho, svd.method = svd.method)
        }
        else {
            out.1 <- Optimize.Theta.multB(Y = Y, maxK = K, B = B, 
                Cov = X, A = A.ine, c = c.ine, D.ker = D.ker, 
                Var.0 = Var.0, tol.rho = tol.rho, max.iter.rho = max.iter.rho, 
                svd.method = svd.method)
        }
        if (!is.null(X) && K > 0) {
            out.1$C <- lapply(out.1$C, function(x, Q.X) {
                if (is.null(x)) {
                  return(NULL)
                }
                return(Q.X \%*\% x)
            }, Q.X = Q.X)
        }
        if (return.all) {
            out$C <- out.1$C
            out$K <- 0:K
            out$rho <- out.1$Rho
            names(out$C) <- out$K
            rownames(out$rho) <- out$K
        }
        else {
            out$C <- out.1$C[[K + 1]]
            out$rho <- out.1$Rho[K + 1, ]
        }
        return(out)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Factor analysis}
