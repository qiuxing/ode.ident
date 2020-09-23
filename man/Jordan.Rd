\name{Jordan}
\alias{Jordan}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generate Jordan blocks with a set of given eigenvalues. }
\description{
  This function creates the Jordan blocks over the field of real numbers
  from a list of eigenvalues.
}
\usage{
Jordan(eigenvals)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eigenvals}{ A list of \code{d} eigenvalues. }
}
\details{
this is a convenient function that creates Jordan blocks (the \code{J}
matrix in a Jordan canonical form) over the field of real numbers from a
list of eigenvalues. Note that for each complex eigenvalue (those with
nonzero imaginary part), this function will automatically add its
complex conjugate to the list of eigenvalues. Therefore, there is no
need to explicitly include both \code{a+b*i} and \code{a-b*i} in the
input. 

The returned value, J matrix, is organized in such way: the first K1
diagonal elements are real eigenvalues; the rest are 2x2 rotational
matrices correspond with pairs of complex eigenvalues.
}
\value{
    \item{J}{The block-diagonal matrix that represents real (1x1 blocks)
      and complex eigenvalues (2x2 blocks) of A. }
}
\references{ \url{https://en.wikipedia.org/wiki/Jordan_normal_form} }
\author{Xing Qiu}
\seealso{\code{\link{eigen}}}
\examples{
ll1 <- c(complex(real=c(-0.1, -0.5), imaginary=3*c(1, 0)))
J1 <- Jordan(ll1)
J1

}                               % end examples.

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{math}
