\name{ptproc.fit}
\alias{ptproc.fit}


\title{
  Fit Point Process Models
}

\description{
  Fit point process models via maximum likelihood.
}

\usage{
ptproc.fit(ppobj, optim.control = list(), method = "Nelder-Mead", 
           alpha = 0, ...) 
}

\arguments{
  \item{ppobj}{A \code{"ptproc"} object containing the model to be fit.}
  \item{optim.control}{An optional named list of tuning parameters to be
    passed to \code{optim}.}
  \item{method}{The optimization method to be used.}
  \item{alpha}{A value for the penalty, if one is needed.}
  \item{...}{Other arguments to be passed to \code{optim}.}
}

\details{
  \code{ptproc.fit} calls the R function \code{optim} to minimize the
  negative log-likelihood.  For minimization, \code{ptproc.fit} uses the
  default method of \code{optim}, which is the Nelder-Mead algorithm.
  This method produces reasonable results for a variety of different
  models.

  Scaling information for the parameters can be passed via the
  \code{parscale} tuning parameter of the \code{optim.control} list.
  This can be very useful for guiding the optimizer to a good solution.

  The value chosen for \code{alpha} will depend on the model and the
  nature of the data.  Generally, one should choose \code{alpha} to be a
  1 or 2 orders of magnitude greater than the typical log-likelihood value.
}

\examples{
data(fires)
ppm <- ptproc(fires[,1], cond.int = hawkes.cond.int,
              params = c(.004, .001, 0),
              condition = penalty(NULL, "any(params < 0)"))
fit <- ptproc.fit(ppm, alpha = 1e+5)
summary(fit)
}






\author{Roger D. Peng}

\keyword{optimize}
