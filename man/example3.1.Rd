\name{example3.1}
\alias{A2}
\alias{tt}
\alias{x0.A}
\alias{x0.B}
\alias{yy1}
\alias{yy2}

\title{Sample solution trajectories from a two-dimensional ODE system. }

\usage{
data(example3.1)
}

\description{
  This dataset contains A sample ODE system (\code{A2}), two initial conditions
  (\code{x0.A} and \code{x0.B}), and the two correponding sets of
  solution trajectories (\code{yy1}, starts with \code{x0.A};
  \code{yy2}, starts with \code{x0.B}). It is called "Example 3.1" in
  our manuscript.
}

\value{
  \item{A2}{ A 2x2 matrix which describes a bivariate linear ODE system.}
  \item{tt}{ The set of 101 time points. }
  \item{x0.A}{ The first initial condition. A2 is practically
  identifiable at \code{x0.A}. }
  \item{x0.B}{ The second initial condition. A2 is not practically
  identifiable at \code{x0.B}. }
  \item{yy1}{ The first set of descrete observations of the solution
  trajectories starting at \code{x0.A}, superimposed with small
  measurement errors. }
  \item{yy2}{ The second set of descrete observations of the solution
  trajectories starting at \code{x0.B}, superimposed with small
  measurement errors. }
}

\references{X. Qiu, T. Xu, B. Soltanalizadeh, and H. Wu. (2020+) Mathematical and Practical Identifiability of High-dimensional Linear Ordinary Differential Equation Systems. Submitted.}
\keyword{datasets}  
