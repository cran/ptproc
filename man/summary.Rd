\name{summary.ptproc}

\alias{summary.ptproc}
\alias{print.summary.ptproc}

\title{
  Summary of point process object
}

\description{
  Print a summary of a \code{"ptproc"} object.
}

\usage{
\method{summary}{ptproc}(object, ...)
\method{print}{summary.ptproc}(x, digits = getOption("digits") - 3, ...)
}

\arguments{
  \item{object,x}{A \code{"ptproc"} object.}
  \item{digits}{The number of digits to be printed for each parameter in
    the summary.}
  \item{...}{Ignored.}
}

\examples{
data(fires)
ppm <- ptproc(fires, cond.int = linear.cond.int,
              params = c(mu = .003, beta1 = .0001,
                         beta2 = .0002, beta3 = .0003))

## Compare against homogeneous Poisson model
summary(ppm)
}

\author{Roger D. Peng}

\keyword{methods}
