\name{twostage1}
\alias{twostage1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{The two-stage estimator of a linear ordinary differential
  equation system based on finite difference method. }
\description{
  This function takes the discrete temporal data as input, and apply
  finite difference method to compute the derivative of the solution
  curves. It then solves an algebraic equation based on the pairwise
  Euclidean inner products to estimate A, the system matrix.
}
\usage{
twostage1(Y, tstep)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{ An \code{dxn}-dimensional matrix of discrete observations. Each
  row is a dimension and each column is a timepoint. }
  \item{tstep}{ The constant difference between two consecutive time
    points (the time step). Currently, this function does not work for
    data with variable time steps. }
}
\details{
  This twostage method is called the simple two-stage method in our
  manuscript.  We also implemented a more accurate two-stage method based on
  smoothing splines as \code{\link{twostage2}} in this R package.
}
\value{
  \item{Ahat}{The estimated system matrix.}
  \item{S}{An \code{nxn}-dimensional matrix that is needed by
    \code{\link{PISAnalysis}}. It represents the pairwise inner product
    between the solution curves and themselves. }
  \item{L}{An \code{nxn}-dimensional matrix that is needed by
    \code{\link{PISAnalysis}}. It represents the pairwise inner product
    between the derivative of solution curves and the solution curves. }
  \item{x0}{The initial condition. }
}
\references{X. Qiu, T. Xu, B. Soltanalizadeh, and H. Wu. (2020+) Identifiability Analysis of Linear Ordinary Differential Equation Systems with a Single Trajectory. Submitted.}
\author{Xing Qiu}
\seealso{\code{\link{fsmooth}}, \code{\link{twostage2}}}
\examples{
## load Example 3.1. In this example, yy1 are discrete and noisy
## observations of (A2, x0.A), which  is practically identifiable.
data("example3.1")
tstep <- tt[2]-tt[1]

myfit1.simple <- twostage1(yy2, tstep)
myfit1.functional <- twostage2(yy1, tt)

## A2 is the true system matrix
round(A2,2)

## simple twostage method is not very accurate
round(myfit1.simple$Ahat,2); round(sum((myfit1.simple$Ahat - A2)^2),2)

## functional twostage method is much more accurate
round(myfit1.functional$Ahat,2); round(sum((myfit1.functional$Ahat - A2)^2),2)


}                               % end examples.

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
