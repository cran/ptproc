\name{penalty}
\alias{penalty}


\title{
  Add penalties to point process object
}

\description{
  Add a penalty to a point process object.
}

\usage{
penalty(code = NULL, condition = FALSE)
}

\arguments{
  \item{code}{A character string representing one or many valid R
    expressions that can be evaluated in an environment containing the
    point process object.}
  \item{condition}{A test to be performed to determine whether or not a
    penalty is incurred.}
}


\details{
  The value of \code{condition} should be a character string
  representing a "test" which will determine whether or not a penalty
  will be returned before evaluation of the log-likelihood.  The penalty
  will be returned if the test evaluates to \code{TRUE}.  The value of
  the penalty is specified in the \code{alpha} parameter to
  \code{ptproc.fit}.

  The \code{code} argument can be used for passing extra code into the
  version of the log-likelihood function passed to the optimizer.
  One can access the full \code{"ptproc"} object via the name `ppobj'.
}


\examples{
data(fires)

## Construct a K = 3 Hawkes cluster model
ppm <- ptproc(fires[,1], cond.int = hawkes.cond.int,
              params = c(.004, .001, rep(0, 3)))

## Restrict parameters to be positive
condition(ppm) <- penalty(code = NULL, condition = "any(params < 0)")

## Use simple linear model
ppm <- ptproc(fires, cond.int = linear.cond.int,
              params = c(.004, .001, 0, 0))

## Evaluate at the "corners" of the domain
extra.code <- paste("ranges <- as.list(as.data.frame(ppobj$ranges))",
                    "corners <- expand.grid(ranges)",
                    "ci <- evalCIF(ppobj, xpts = corners)", sep = ";")

## Add penalty
condition(ppm) <- penalty(code = extra.code, condition = "any(ci < 0)")
}

\author{Roger D. Peng}

\keyword{utilities}
