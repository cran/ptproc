\name{ptproc.sim}
\alias{ptproc.sim}

\title{
  Simulate a point process
}

\description{
  Simulate a point process, given a model.
}

\usage{
ptproc.sim(ppobj, M)
}

\arguments{
  \item{ppobj}{A \code{"ptproc"} object containing the conditional
    intensity model to be simulated.}
  \item{M}{The maximum of the given conditional intensity over the
    range.}
}


\details{
  The value of \code{M} should be such that M >= lambda(x) for all x in
  the domain (where lambda(x) is the conditional intensity function).
}


\examples{
## Simulate from a linear model
ppsim <- ptproc(pts = NA, cond.int = linear.cond.int,
                params = c(100, 10, -50, 50), ndim = 3,
                ranges = cbind(c(0,2),c(0,2),c(0,2)),
                is.pois = TRUE)
sim <- ptproc.sim(ppsim, M = 330)
colnames(sim) <- c("Time", "X", "Y")
pairs(sim, pch=".")
}


\author{Roger D. Peng}
\keyword{models}
