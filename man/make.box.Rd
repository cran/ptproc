\name{make.box}

\alias{make.box}

\title{
  Get coordinates of a bounding box
}

\description{
  Return a 2-D polygon representing the box specified by the
  \code{ranges} element of a \code{"ptproc"} object.  This can be useful
  when applying spatial statistical tests (such as the K function).
}

\usage{
make.box(ppobj, idx = 1:2)
}

\arguments{
  \item{ppobj}{A \code{"ptproc"} object.}
  \item{idx}{A vector of length 2 specifying the dimensions of the data
    for which the box should be constructed.}
}

\author{Roger D. Peng}

\keyword{utilities}
