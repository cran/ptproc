\name{CIF.util}

\alias{CIF.util}
\alias{evalCIF}
\alias{integrateCIF}

\title{
  Functions for accessing the conditional intensity
}

\description{
  Functions for evaluating and integrating the conditional intensity.
}

\usage{
integrateCIF(ppobj, TT = ppobj$ranges, ...)
evalCIF(ppobj, xpts = ppobj$pts, ...)
}

\arguments{
  \item{ppobj}{A \code{"ptproc"} object.}
  \item{TT}{A matrix containing the ranges of integration in each
    dimension.  If not specified, the \code{ranges} element of
    \code{ppobj} is used.}
  \item{xpts}{A matrix of points where the conditional intensity will be
    evaluated.}
  \item{...}{Other arguments passed to the conditional intensity
    function.}
}

\details{
  These functions simply call the evaluation and integration portions of
  the (possibly user-defined) conditional intensity function.
  
  \code{evalCIF} evaluates the conditional intensity specified with the
  model at a matrix of points.  By default, it evaluates at the data
  points.  \code{integrateCIF} returns the integral
  of the conditional intensity over the given ranges in \code{TT}, the
  default being the entire domain.  
}
  
\author{Roger D. Peng}

\keyword{utilities}
