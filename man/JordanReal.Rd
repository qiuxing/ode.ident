\name{JordanReal}
\alias{JordanReal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Jordan canonical decomposition of a real-valued matrix over the
  field of real numbers. }
\description{
  This function computes the Jordan canonical form of a a real-valued
  matrix A over the field of real numbers. 
}
\usage{
JordanReal(A,tol=1e-6)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{ A \code{dxd}-dimensional real-valued matrix. }
  \item{tol}{ The numerical tolerance level. Due to inevitable 
    limitations of numerical computations, it is possible that all
    eigenvalues computed by \code{eigen(A)} have some small imaginary
    component even if some of them are clearly real. In
    \code{JordanReal}, if the absolute value of the imaginary part of an eigenvalue is less than \code{tol}, we treat it as a real eigenvalue. }
}
\details{
this is a convenient function that computes the Jordan canonical form
(JCF) of a real-valued matrix A over the field of real numbers. Please
note that this function is \strong{not} numerically stable for all
matrices (esp. large matrices). It assumes that there is no nilpotent
cells in the J.C.F.  Returned values: K1 is the number of real eigenvalues; K2 is the number of
pairs of complex eigenvalues. J is the pxp semi-diagonal matrix of
Jordan blocks; Qmat is the matrix of generalized eigenvectors. Note
that J is always organized in such way: the first K1 diagonal
elements are real eigenvalues; the rest are 2x2 rotational matrices
correspond with pairs of complex eigenvalues.
}
\value{
    \item{J}{The block-diagonal matrix that represents real (1x1 blocks)
      and complex eigenvalues (2x2 blocks) of A. }
    \item{Qmat}{The matrix of eigenvectors and semi-eigenvectors (for
      complex eigenvalues). }
    \item{Qinv}{The inverse of \code{Qmat}.}
    \item{K1}{Number of real eigenvalues.}
    \item{K2}{Number of pairs of complex eigenvalues. So \code{d=K1+2*K2}. }

}
\references{ \url{https://en.wikipedia.org/wiki/Jordan_normal_form} }
\author{Xing Qiu}
\seealso{\code{\link{eigen}}}
\examples{
A1  <- matrix(c(0, 1, -1,
                2, 0, 0,
                3, 1, 0), 3, byrow=TRUE)
jcf1 <- JordanReal(A1)
jcf1
}                               % end examples.

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{math}
