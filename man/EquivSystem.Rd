\name{EquivSystem}
\alias{EquivSystem}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Equivalent system matrices for an unidentifiable linear ODE system. }
\description{
  This function generates an equivalent system matrix for 
  a linear ODE system that is unidentifiable at \code{x0} based on a
  given coefficient matrix \code{D}.
}
\usage{EquivSystem(A, x0, Dmat, ...)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{ The original system matrix. }
  \item{x0}{ The initial condition.}
  \item{Dmat}{ When A is not identifiable at \code{x0}, the class of all
    equivalent system matrices is an affine space with nonzero degrees
    of freedom. Dmat is the coefficient matrix that determines the value
    of one specific equivalent system matrix in that affine space.}
  \item{...}{ Additional parameters used by function
    \code{ICISAnalysis(A, x0, ...)}.}
}
\details{ EquivSystem computes \code{[A, x0]}, the equivalence class of
A given x0, and then creates an equivalent system A (for the same x0)
based on the specified difference matrix Dmat (with the same
dimensionality of A).  Remarks: a) when Dmat is the zero matrix,
EquivSystem(A, x0, 0)=A; b) when A is *identifiable* at x0,
EquivSystem(A, x0, Dmat)=A for all input Dmat.  For unidentifiable
cases, using a random Dmat a.s. will lead to a different system Atilde
that is equivalent to A. }
\value{
  \item{Atilde}{A system matrix that is equivalent to A at \code{x0}.}
}
\references{ X. Qiu, T. Xu, B. Soltanalizadeh, and H. Wu. (2020+)
  Mathematical and Practical Identifiability of High-dimensional Linear
  Ordinary Differential Equation Systems. Submitted.}
\author{Xing Qiu}
\seealso{ \code{\link{ICISAnalysis}}, \code{\link{PISAnalysis}} }
\examples{

A1  <- matrix(c(0, 1, -1,
                2, 0, 0,
                3, 1, 0), 3, byrow=TRUE)
jcf1 <- JordanReal(A1); J1 <- jcf1$J; Qmat1 <- jcf1$Qmat; Qinv1 <- jcf1$Qinv

## A1 is identifiable at x0.a but not identifiable at x0.b
x0.a <- Qmat1 \%*\% c(2, -1, 0)
x0.b <- Qmat1 \%*\% c(0, -2, 3)

## Using a random D matrix to create an equivalent system matrix
D1 <- matrix(rnorm(3*3), 3)

## Note that, because A1 is identifiable at x0.a, adding a random D1
## matrix will not change A1 at all.
B1 <- EquivSystem(A1, x0.a, Dmat=D1)
all.equal(A1, B1) #yes

## Now try x0.b. This time, we will get a system matrix B2 that is
## totally different from A1.
B2 <- EquivSystem(A1, x0.b, Dmat=D1)
B2; A1

}                               % end examples.

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{math}
