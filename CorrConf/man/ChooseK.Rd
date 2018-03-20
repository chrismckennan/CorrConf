\name{ChooseK}
\alias{ChooseK}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate the latent factor dimension K
}
\description{
Estimate the latent factor dimension K, the number of columns in the unobserved covariate matrix, C. The user should supply the expression matrix, Y, the observed covariates to be included in the model, Cov, the maximum number of K's to consider, maxK, and the list (or matrix if b = 1) of relatedness matrices. If the user does not supply the identity I_n or if I_n is not in the span of the other relatedness matrices, I_n will be automatically added to the list of relatedness matrices.
}
\usage{
ChooseK(Y, Cov = NULL, maxK = 20, B = NULL, nFolds = 10, simpleDelta = F, A.ine = NULL, c.ine = NULL, A.equ = NULL, Var.0 = NULL, tol.rho = 0.001, max.iter.rho = 15, svd.method = "fast", plotit = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
A p x n data matrix; p = number of units (i.e. genes), n = sample size
}
  \item{Cov}{
A n x (d+r) covariate matrix, including the d covariates of interest and r nuisance covaraites (like the intercept)
}
  \item{maxK}{
Maximum latent dimension K to consider
}
  \item{B}{
A list of length b of n x n positive semi-definite relatedness matrices if b > 1 or a matrix if b = 1. If the user does not supply the identity I_n or if I_n is not in the span of the other relatedness matrices, I_n will be automatically added to the list of relatedness matrices. If left unspecified, it is assumed all individuals (i.e. samples) are uncorrelated.
}
  \item{nFolds}{
Number of folds in cross validation.
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
A #equality-constraints x b matrix giving the equality constraints on the variance multipliers. If left unspecified, it defaults to assuming there are no equality constraints. If this is specified, the user MUST define \code{Var.0} and it must be in the feasible region.
}
  \item{Var.0}{
A starting point for tau. If there are equality constraints, the user MUST specifiy a starting point within the feasible region.
}
  \item{tol.rho}{
Each iteration of sequential PCA completes when || tau_j/||tau_j||_2 - tau_\{j-1\}/||tau_\{j-2\}||_2 ||_2 <= b * tol.rho
}
  \item{max.iter.rho}{
Maximum number of iterations j for each sequential PCA iteration
}
  \item{svd.method}{
"fast": uses the package \code{irlba} for svd; any other options uses the R default.
}
  \item{plotit}{
If TRUE, it plots the leave-one-out cross-validated error as a function of K
}
}
\details{
The true K may not be the one with the lowest cross-validation error. It is recommended the user try different values around the K with the minimum LOO-XV error to ensure their results are not too dependent on the choice of K.
}
\value{
A data frame
\item{K}{A vector of latent dimensions considered}
\item{LOO.XV}{A vector of the average leave-one-out cross validation for each \code{K} considered}
\item{K.hat}{The \code{K} that gives the minimum leave-one-out cross validation}
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
function (Y, Cov = NULL, maxK = 20, B = NULL, nFolds = 10, simpleDelta = T, 
    A.ine = NULL, c.ine = NULL, A.equ = NULL, Var.0 = NULL, tol.rho = 0.001, 
    max.iter.rho = 15, svd.method = "fast", plotit = T) 
{
    if (is.null(B)) {
        return(ChooseK_NoB(Y = Y, X = Cov, maxK = maxK, nFolds = nFolds, 
            simpleDelta = simpleDelta, max.iter.svd = 3, svd.method = svd.method, 
            plotit = plotit))
    }
    if (is.matrix(B) || (is.list(B) && length(B) == 1)) {
        if (is.list(B)) {
            B <- B[[1]]
        }
        if (simpleDelta) {
            return(ChooseK_parallel.simrho(Y = Y, X = Cov, maxK = maxK, 
                B = B, nFolds = nFolds, tol.rho = tol.rho, max.iter.rho = max.iter.rho, 
                svd.method = svd.method, plotit = plotit))
        }
        else {
            return(ChooseK_parallel(Y = Y, X = Cov, maxK = maxK, 
                B = B, nFolds = nFolds, tol.rho = tol.rho, max.iter.rho = max.iter.rho, 
                svd.method = svd.method, plotit = plotit))
        }
    }
    if (is.list(B) && length(B) > 1) {
        B <- IncludeIdent(B)
        D.ker <- CreateD.ker(A.equ)
        if (simpleDelta) {
            return(ChooseK_parallel.multB.simrho(Y = Y, X = Cov, 
                maxK = maxK, B = B, nFolds = nFolds, A.lin = A.ine, 
                c.lin = c.ine, D.ker = D.ker, Var.0 = Var.0, 
                tol.rho = tol.rho, max.iter.rho = max.iter.rho, 
                svd.method = svd.method, plotit = plotit))
        }
        else {
            return(ChooseK_parallel.multB(Y = Y, X = Cov, maxK = maxK, 
                B = B, nFolds = nFolds, A.lin = A.ine, c.lin = c.ine, 
                D.ker = D.ker, Var.0 = Var.0, tol.rho = tol.rho, 
                max.iter.rho = max.iter.rho, svd.method = svd.method, 
                plotit = plotit))
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Leave-one-out cross validation}
