\name{twostage2}
\alias{twostage2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{The two-stage estimator of a linear ordinary differential
  equation system based on smoothing splines. }
\description{
  This function takes the discrete temporal data as input, and first
  apply roughness penalized smoothing spline to represent these data as
  smooth functions (basis splines). It then solves an algebraic equation
  based on the pairwise inner products to estimate A, the system matrix.
}
\usage{
twostage2(Y, Ts, nord=4, rough.pen=1e-3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{ An \code{dxn}-dimensional matrix of discrete observations. Each
  row is a dimension and each column is a timepoint. }
  \item{Ts}{ An \code{n}-dimensional vector of time points in ascending
  order. Our recommendation is to standardize those time points so that
  \code{Ts[1]}=0 and \code{Ts[length(Ts)]}=1. }
  \item{nord}{ The order of b-splines, which is one higher than their
  degree. The default of 4 gives cubic splines. }
  \item{rough.pen}{Roughness penalty used in smoothing spline. Its
  default value is 0.001, which is a reasonable value if: (a) the range
  of\code{Ts} is relatively small (e.g., from 0 to 1 as we recommended), and
  (b) the signal-to-noise level of the data is relatively small, i.e.,
  it is easy to see the overall temporal trend from the discrete data by
  visual examination. }
}
\details{
  This twostage method is called the functional two-stage method in our
  manuscript.  We also implemented a simpler two-stage method based on
  finite difference as \code{\link{twostage1}} in this R package.
}
\value{
  \item{Ahat}{The estimated system matrix.}
  \item{S}{An \code{nxn}-dimensional matrix that is needed by
    \code{\link{PISAnalysis}}. It represents the pairwise inner product
    between the solution curves and themselves. }
  \item{L}{An \code{nxn}-dimensional matrix that is needed by
    \code{\link{PISAnalysis}}. It represents the pairwise inner product
    between the derivative of solution curves and the solution curves. }
  \item{xt.hat}{The smoothed curves used by the functional two-stage
    method. }
  \item{x0}{The estimated initial condition. }
}
\references{X. Qiu, T. Xu, B. Soltanalizadeh, and H. Wu. (2020+) Identifiability Analysis of Linear Ordinary Differential Equation Systems with a Single Trajectory. Submitted.}
\author{Xing Qiu}
\seealso{\code{\link{fsmooth}}, \code{\link{twostage1}}}
\examples{
## load Example 3.1. In this example, yy1 are discrete and noisy
## observations of (A2, x0.A), and yy2 are observations of (A2,x0.B).
## The first case is practically identifiable but the second is not.
data("example3.1")

## To demonstrate that PIS is useful, we will use the functional
## two-stage method to estimate A2
myfit1 <- twostage2(yy1, tt)
S <- myfit1$S; L <- myfit1$L
myfit2 <- twostage2(yy2, tt)

## PISs of the first case are all relatively large
## (good practical identifiability) 
PISs1 <- PISAnalysis(yy1, S, L); PISs1
## PISs of the second case are all relatively small
## (bad practical identifiability) 
PISs2 <- PISAnalysis(yy2, S, L); PISs2

## A2 is the true system matrix
round(A2,2)
## The first case is a good estimate of A2
round(myfit1$Ahat,2); round(sum((myfit1$Ahat - A2)^2),2)
## The second case is a bad estimate of A2 due to
## identifiability issues
round(myfit2$Ahat,2); round(sum((myfit2$Ahat - A2)^2),2)

}                               % end examples.

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
