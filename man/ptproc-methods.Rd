\name{ptproc.methods}

\alias{ptproc.methods}

\alias{logLik.ptproc}
\alias{print.ptproc}
\alias{residuals.ptproc}

\title{
  Methods for point process objects
}

\description{
  Some methods for handling \code{"ptproc"} objects.
}

\usage{
\method{logLik}{ptproc}(object, negative = FALSE, ...)
\method{print}{ptproc}(x, digits = getOption("digits") - 3, ...)
\method{residuals}{ptproc}(object, type = c("ordinary", "approx"),
                           m = NULL, K = NULL, ...)
}

\arguments{
  \item{object, x}{An object of class \code{"ptproc"}.}
  \item{negative}{If set to \code{TRUE} the negative log-likelihood is
    returned.}
  \item{digits}{The number of digits to be printed for the parameters in
    the \code{"ptproc"} object.}
  \item{type}{The type of residual process to produce.}
  \item{m}{The minimum value of the conditional intensity over the domain.}
  \item{K}{The subsample size to be drawn when approximate residuals are
  chosen.}
  \item{...}{Other arguments passed to and from other methods.}
}


\examples{
## None right now

data(fires)
ppm <- ptproc(fires[,1], cond.int = hawkes.cond.int,
              params = c(mu = .004, C = .001, a = 0.1),
              condition = penalty(NULL, "any(params < 0)"))
print(ppm)
logLik(ppm)
}

\author{Roger D. Peng}


\keyword{methods}
