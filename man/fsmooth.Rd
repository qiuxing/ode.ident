\name{fsmooth}
\alias{fsmooth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A wrapper for smoothing splines}
\description{
  This function is a wrapper to smooth a set of multidimensional discrete timecourse data.

  timecourse data; Ts are the time points; rough.pen is the roughness
  penalty used in smoothing.  Note that we must be very careful in
  selecting the "right" penalty for the real data.  n.plot: number of
  functions to be plotted in the figure (having too many curves in
  one figure is not a good idea).
}
\usage{
fsmooth(Y, Ts, rough.pen=.001, norder=4, plot=FALSE, n.plot=5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{ An \code{dxn}-dimensional matrix of discrete observations. Each
  row is a dimension and each column is a timepoint. }
  \item{Ts}{ An \code{n}-dimensional vector of time points in ascending
  order. Our recommendation is to standardize those time points so that
  \code{Ts[1]}=0 and \code{Ts[length(Ts)]}=1. }
  \item{rough.pen}{Roughness penalty used in smoothing spline. Its
  default value is 0.001, which is a reasonable value if: (a) the range
  of\code{Ts} is relatively small (e.g., from 0 to 1 as we recommended), and
  (b) the signal-to-noise level of the data is relatively small, i.e.,
  it is easy to see the overall temporal trend from the discrete data by
  visual examination. }
  \item{norder}{ The order of b-splines, which is one higher than their
  degree. The default of 4 gives cubic splines. }
  \item{plot}{If \code{plot=TRUE}, A plot of the discrete data and 
  smoothed curves in the first \code{n.plot} dimensions will be plotted. The default value is \code{FALSE}. }
  \item{n.plot}{ Number of dimensions to be used in the plot. }
}
\details{
  This function is a wrapper to smooth a set of multidimensional
  discrete timecourse data. It calls functions
  \code{create.bspline.basis()}, \code{fdPar()}, and
  \code{smooth.basis()} in R package \code{fda} to perform the
  computation.  
}
\value{
  \item{Xt}{A list of smoothed cubic splines.}
}
\references{ Ramsay, James O., Hooker, Giles, and Graves, Spencer (2009), Functional data analysis with R and Matlab, Springer, New York.}
\author{Xing Qiu}
\seealso{\code{\link{create.bspline.basis}}, \code{\link{fdPar}}, \code{\link{smooth.basis}}}
\examples{

## load Example 3.1
data("example3.1")

## smoothing
yt1 <- fsmooth(yy1, tt, rough.pen=1e-3, plot=FALSE)
yt2 <- fsmooth(yy2, tt, rough.pen=1e-3, plot=TRUE)

}                               % end examples.

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{smooth}
