\name{ptproc}
\alias{ptproc}

\title{
  Point process constructor
}

\description{
  Specify and create a \code{ptproc} object.
}

\usage{
ptproc(pts = NA, cond.int, params, fixed.params = rep(NA, length(params)),
       condition = NULL, ranges = NULL, data = NULL, ndim = NULL,
       is.pois = FALSE)
}

\arguments{
  \item{pts}{A n x p matrix of data points where n is the number of
    points and p is the dimension of the data.}
  \item{cond.int}{A conditional intensity function.}
  \item{params}{A vector of parameters for the point process model.}
  \item{fixed.params}{A vector of \code{NA} and non-\code{NA} elements
    indicating which parameters are fixed.}
  \item{condition}{A valid unevaluated R expression (see Details).}
  \item{ranges}{A 2 x p matrix indicating the ranges for each dimension
    in the domain.}
  \item{data}{Optional extra information for the point process model.}
  \item{ndim}{The dimension of the data.}
  \item{is.pois}{Is the point process model a Poisson model?}
}

\details{
  For simulating a point process, \code{pts} can be specified as
  \code{NA} since no data points are needed for simulation.  
  
  The \code{condition} element is used for including penalties in the
  model for the evaluation of the log-likelihood.  See the help page for
  \code{penalty} for more information.

  At this point the domain of the point process can only be a "box", the
  boundaries of which are specified in \code{ranges}.
  The (1, j) element of the \code{ranges} element contains the lower
  bound of the range in the jth dimension.  Similarly, the (2, j)
  element of \code{ranges} contains the upper bound of the range in the
  jth dimension.

  The \code{data} element can be used to include other information into
  the point process model.  Often this will be in the form of a (named)
  list containing information about covariates, marks, or perhaps
  preprocessed values.

  The \code{is.pois} element does not need to be set and defaults to
  \code{FALSE}.  However, if the model is indeed Poisson and
  \code{is.pois} is set to \code{TRUE}, then this may speed up
  simulation of the model.
}

\seealso{
  The \code{"ptproc"} object has \code{\link{print}} and
  \code{\link{summary}} methods.

  For example conditional intensity functions see
  \code{\link{cond.int}}.
}


\author{Roger D. Peng}

\examples{
data(fires)

## Use a linear conditional intensity model
pp.model <- ptproc(pts = fires, cond.int = linear.cond.int,
                   params = c(.004, 0, 0, 0))

## Use a Hawkes model and keep a parameter fixed
pp.model <- ptproc(pts = fires[,1], cond.int = hawkes.cond.int,
                   params = c(.004, .004, rep(0, 2)),
                   fixed.params = c(NA, 1, NA, NA))
}

\keyword{models}
