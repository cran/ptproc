\name{cond.int}
\alias{cond.int}
\alias{linear.cond.int}
\alias{hawkes.cond.int}
\alias{hPois.cond.int}

\title{
  Conditional intensity functions
}

\description{
  Example point process conditional intensity functions.
}

\usage{
linear.cond.int(params, eval.pts, pts = NA, data = NULL, TT = NULL)
hawkes.cond.int(params, eval.pts, pts = NA, data = NULL, TT = NULL)
hPois.cond.int(params, eval.pts, pts = NA, data = NULL, TT = NULL)
}

\arguments{
  \item{params}{A vector of parameters for a point process model.}
  \item{eval.pts}{A matrix of points at which the conditional intensity
    should be evaluated.}
  \item{pts}{A matrix containing original data.}
  \item{data}{Other information that may be needed to evaluate or
    integrate the conditional intensity.}
  \item{TT}{A matrix containing the ranges of integration for each
    dimension.}
}


\details{
  \code{linear.cond.int} implements a simple space-time linear model
  with a background parameter mu and a parameter for the time and
  spatial dimensions.

  \code{hawkes.cond.int} implements a one-dimensional Hawkes-type
  cluster model with a Laguerre-type polynomial as the trigger
  function.  See the reference for more details on this model.

  \code{hPois.cond.int} is a simple homogeneous Poisson model.  It is
  unlikely to be useful but is included as yet another example.
}

\references{
  Ogata, Y. (1988) Statistical models for earthquake occurrences and
  residual analysis for point processes. \emph{JASA}, \bold{83}, 9--27.
}
   

\examples{
## See ?ptproc for examples
}



\author{Roger D. Peng}

\keyword{models}
