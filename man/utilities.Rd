\name{ptproc.util}

\alias{ptproc.util}
\alias{fixed}
\alias{fixed<-}
\alias{params}
\alias{params<-}
\alias{condition}
\alias{condition<-}



\title{
  Utilities functions for ptproc objects
}

\description{
  Utility functions for accessing/modifying \code{"ptproc"} objects.
}

\usage{
fixed(x)
fixed(x) <- value
params(x)
params(x) <- value
condition(x)
condition(x) <- value
}

\arguments{
  \item{x}{A \code{"ptproc"} object.}
  \item{value}{For \code{fixed} and \code{params}, a vector of numbers.
    For \code{condition}, a valid R expression object.}
}

\details{
  These functions can be used to access and modify various elements of
  the \code{"ptproc"} object.
}

\examples{
data(fires)
ppm <- ptproc(pts = fires, cond.int = linear.cond.int,
              params = c(.004, 0, 0, 0))
## Fix the background
fixed(ppm)[1] <- 1

## Change the y-direction parameter
params(ppm)[4] <- .001

## Add a condition to restrict parameters positive
condition(ppm) <- penalty(code = NULL, condition = "any(params < 0)")
}




\keyword{attribute}
