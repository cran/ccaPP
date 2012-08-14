\name{fastMedian}
\alias{fastMedian}
\title{Fast implementation of the median}
\usage{
  fastMedian(x)
}
\arguments{
  \item{x}{a numeric vector.}
}
\value{
  The sample median.
}
\description{
  Compute the sample median with a fast C++ implementation.
}
\note{
  Functionality for removing observations with missing
  values is currently not implemented.
}
\examples{
set.seed(1234)  # for reproducibility
x <- rnorm(100)
fastMedian(x)
}
\author{
  Andreas Alfons
}
\seealso{
  \code{\link{fastMAD}}, \code{\link[stats]{median}}
}
\keyword{multivariate}
\keyword{robust}

