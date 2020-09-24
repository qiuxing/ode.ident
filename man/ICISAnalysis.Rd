\name{ICISAnalysis}
\alias{ICISAnalysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Perform identifiability analysis of an ODE system A at the
  initial value \code{x0}.}
\description{
  this function computes the initial condition-based identifiability
  score (ICIS) and a few related quantitative measures of
  identifiability of an ODE system at initial conditon x0.   
}
\usage{
ICISAnalysis(A, x0, n.digits=3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{ The original system matrix. }
  \item{x0}{ The initial condition.}
  \item{n.digits}{ Number of significant digits used in computing
    \code{w0k.norm} hence \code{ICIS}. }
}
\details{
  The initial condition-based identifiability score ICIS is a
  quantitative measure of identifiability of an ODE system A at the
  initial value \code{x0}. If ICIS equals zero (the significant digits
  of ICIS is less than \code{n.digits}), A is not mathematically
  identifiable at \code{x0}. Besides, A may be unidentifiable due to
  having repeated eigenvalues, which will also be checked by \code{ICISAnalysis}.
}
\value{
  \item{Identifiable}{A logic variable, TRUE means identifiable, FALSE
    means unidentifiable.}
  \item{Ident1}{A logic variable, TRUE means ICIS is greater than zero,
    FALSE means ICIS is zero (up to \code{n.digits} of significant
    digits) therefore A is not identifiable at \code{x0}.}
  \item{Ident2}{A logic variable, TRUE means there is no repeated
    eigenvalue of A (up to \code{n.digits} of significant
    digits), FALSE means there are repeated eigenvalues of A, therefore
    A is not identifiable (for any \code{x0}). }
  \item{ICIS}{The initial condition-based identifiability score. Larger
    value of ICIS implies better practical identifiability.  }
  \item{w0k.norm}{ A vector of K elements. \code{w0k.norm[k]} represents the
    projection of \code{x0} on the linear space spanned by the kth
    \emph{invariant subspace} of A. Note that ICIS is defined as the
    minimum of \code{w0k.norm}. }
  \item{I0}{A \code{dxd}-dimensional diagonal matrix with either zero or
    one on the diagonal. If the ith diagonal entry is one, A is
    not identifiable on the invariant subspace corresponds with the ith
    eigenvector. }
  \item{Iplus}{ A \code{dxd}-dimensional diagonal matrix that is the
    complement of \code{I0}. \eqn{I_{+} := I_{d\times d} - I_{0}}. }
  \item{Jordan}{ Some useful information provided by the Jordan
    decomposition of A on the field of real numbers.
    \describe{
      \item{J}{The block-diagonal matrix that represents real (1x1 blocks)
        and complex eigenvalues (2x2 blocks) of A. }
      \item{Qmat}{The matrix of eigenvectors and semi-eigenvectors (for
        complex eigenvalues). }
      \item{Qinv}{The inverse of \code{Qmat}.}
      \item{K1}{Number of real eigenvalues.}
      \item{K2}{Number of pairs of complex eigenvalues. So \code{d=K1+2*K2}. }
      \item{Lgap}{The smallest gap between eigenvalues. If this quantity
        is zero (up to \code{n.digits} of significant digits), A is not
        identifiable for all \code{x0}. }
      \item{RepEigenIdx}{ A list of the indices of eigenvalues that have
        repeated values. }
    }
  }
}
\references{ X. Qiu, T. Xu, B. Soltanalizadeh, and H. Wu. (2020+) Identifiability Analysis of Linear Ordinary Differential Equation Systems with a Single Trajectory. Submitted. }
\author{Xing Qiu}
\seealso{\code{\link{PISAnalysis}}}
\examples{

A1  <- matrix(c(0, 1, -1,
                2, 0, 0,
                3, 1, 0), 3, byrow=TRUE)
jcf1 <- JordanReal(A1); J1 <- jcf1$J; Qmat1 <- jcf1$Qmat; Qinv1 <- jcf1$Qinv

## A1 is identifiable at x0.a but not identifiable at x0.b
x0.a <- Qmat1 \%*\% c(2, -1, 0)
x0.b <- Qmat1 \%*\% c(0, -2, 3)

## now check the identifiability of A1 at x0.a and x0.b
ICISAnalysis(A1, x0.a) #yes
ICISAnalysis(A1, x0.b) #no because ICIS==0; I0 is not a zero matrix

}                               % end examples.

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{math}
