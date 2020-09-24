\name{PISAnalysis}
\alias{PISAnalysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Perform identifiability analysis of an ODE system A with the
  observed discrete data matrix.}
\description{
  This function computes three types of practical identifiability
  scores, the smoothed condition number (SCN), the practical
  identifiability score (PIS), and Stanhope's kappa statistic. 
}
\usage{
PISAnalysis(Y,S,L)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{ An \code{dxn}-dimensional matrix of discrete observations. Each
  row is a dimension and each column is a timepoint. }
  \item{S}{ See details.}
  \item{L}{ See details.}
}
\details{
  This function computes three types of identifiability scores, SCN,
  PIS, and Stanhope's kappa. Y is the matrix of discrete observations,
  likely with measurement error. S, L are two \code{nxn}-dimensional
  matrices that represents the estimated inner product matrices between the
  solution curves, and between the derivatives of solution curves and
  the solution curves. In other words

  \deqn{ \hat{\Sigma}_{\mathbf{x}\mathbf{x}} = Y S Y', \quad \hat{\Sigma}_{D\mathbf{x}, \mathbf{x}} = Y L Y'. }

  These two matrices are provided by function \code{twostage2}. 

}
\value{
  \item{SCN}{The smoothed condition number (SCN).}
  \item{PIS}{The practical identifiability score (PIS).}
  \item{kappa}{Stanhope's kappa statistic.}
}
\references{ \enumerate{
    \item Stanhope, S., Rubin, J. E., & Swigon,
    D. (2014). Identifiability of linear and linear-in-parameters
    dynamical systems from a single trajectory. SIAM Journal on Applied
    Dynamical Systems, 13(4), 1792-1815.
  \item X. Qiu, T. Xu, B. Soltanalizadeh, and H. Wu. (2020+)
Identifiability Analysis of Linear Ordinary Differential Equation Systems with a Single Trajectory. Submitted.
}}
\author{Xing Qiu}
\seealso{\code{\link{ICISAnalysis}}, \code{\link{twostage2}}}
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
\keyword{math}
