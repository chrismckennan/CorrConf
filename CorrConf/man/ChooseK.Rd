\name{ChooseK}
\alias{ChooseK}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate the latent factor dimension K
}
\description{
Estimate the latent factor dimension K, i.e. the number of columns in the unobserved covariate matrix, C. The user should supply the expression matrix, Y, the observed covariates to be included in the model, Cov, the maximum number of K's to consider, maxK, and the list (or matrix if b = 1) of relatedness matrices.
}
\usage{
ChooseK(Y, Cov = NULL, maxK = 20, B = NULL, nFolds = 10, simpleDelta = F, tol.rho = 0.001, max.iter.rho = 15, svd.method = "fast")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
A p x n data matrix; p = number of units (i.e. genes), n = sample size
}
  \item{Cov}{
A n x (d+r) covariate matrix, including the d covariates of interest and r nuisance covaraites
}
  \item{maxK}{
Maximum latent dimension K to consider
}
  \item{B}{
A list of length b of n x n positive semi-definite relatedness matrices if b > 1 or a matrix if b = 1. If unspecified, it is assumed all individuals (i.e. samples) are uncorrelated.
}
  \item{nFolds}{
Number of folds in cross validation.
}
  \item{simpleDelta}{
If \code{TRUE}, then it assumes Delta = delta^2 * I_p. Otherwise, Delta = diag(delta_1^2,...,delta_p^2). We recommend the user assign this \code{TRUE} for datasets with large p and n, i.e. 800K methylation array data. Otherwise, the program will be prohibitively slow.
}
  \item{tol.rho}{
Each iteration of sequential PCA completes when ||rho_j - rho_{j-1}||_2 <= b * tol.rho
}
  \item{max.iter.rho}{
Maximum number of iterations j for each sequential PCA iteration
}
  \item{svd.method}{
"fast": uses the package \code{irlba} for svd; any other options uses the R default.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A data frame
\item{K}{A vector of latent dimensions considered}
\item{LOO.XV}{A vector of the average leave-one-out cross validation for each \code{K} considered}
\item{K.hat}{The \code{K} that gives the minimum leave-one-out cross validation}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Chris McKennan
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (Y, Cov = NULL, maxK = 20, B = NULL, nFolds = 10, simpleDelta = F, 
    tol.rho = 0.001, max.iter.rho = 15, svd.method = "fast") 
{
    if (is.null(B)) {
        return(ChooseK_NoB(Y = Y, X = Cov, maxK = maxK, nFolds = nFolds, 
            max.iter.svd = 3, svd.method = "fast"))
    }
    if (is.matrix(B) || (is.list(B) && length(B) == 1)) {
        if (is.list(B)) {
            B <- B[[1]]
        }
        if (simpleDelta) {
            return(ChooseK_parallel.simrho(Y = Y, X = Cov, maxK = maxK, 
                B = B, nFolds = nFolds, tol.rho = tol.rho, max.iter.rho = max.iter.rho, 
                svd.method = svd.method))
        }
        else {
            return(ChooseK_parallel(Y = Y, X = Cov, maxK = maxK, 
                B = B, nFolds = nFolds, tol.rho = tol.rho, max.iter.rho = max.iter.rho, 
                svd.method = svd.method))
        }
    }
    if (is.list(B) && length(B) > 1) {
        if (simpleDelta) {
            return(ChooseK_parallel.multB.simrho(Y = Y, X = Cov, 
                maxK = maxK, B = B, nFolds = nFolds, tol.rho = tol.rho, 
                max.iter.rho = max.iter.rho, svd.method = svd.method))
        }
        else {
            return(ChooseK_parallel.multB(Y = Y, X = Cov, maxK = maxK, 
                B = B, nFolds = nFolds, tol.rho = tol.rho, max.iter.rho = max.iter.rho, 
                svd.method = svd.method))
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Leave-one-out cross validation}
