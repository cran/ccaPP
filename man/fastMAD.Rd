\name{fastMAD}
\alias{fastMAD}
\title{Fast implementation of the median absolute deviation}
\usage{
  fastMAD(x, constant = 1.4826)
}
\arguments{
  \item{x}{a numeric vector.}

  \item{constant}{a numeric multiplication factor.  The
  default value yields consistency at the normal model.}
}
\value{
  A list with the following components:

  \item{center}{a numeric value giving the sample median.}

  \item{MAD}{a numeric value giving the median absolute
  deviation.}
}
\description{
  Compute the median absolute deviation with a fast C++
  implementation.  By default, a multiplication factor is
  applied for consistency at the normal model.
}
\note{
  Functionality for removing observations with missing
  values is currently not implemented.
}
\examples{
set.seed(1234)  # for reproducibility
x <- rnorm(100)
fastMAD(x)
}
\author{
  Andreas Alfons
}
\seealso{
  \code{\link{fastMedian}}, \code{\link[stats]{mad}}
}
\keyword{multivariate}
\keyword{robust}

