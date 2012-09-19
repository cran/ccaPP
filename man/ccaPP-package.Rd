\name{ccaPP-package}
\alias{ccaPP-package}
\alias{ccaPP}
\docType{package}
\title{
(Robust) canonical correlation analysis via projection pursuit
}
\description{
Canonical correlation analysis and maximum correlation via 
    projection pursuit, as well as fast implementations of correlation
    estimators, with a focus on robust and nonparametric methods.
}
\details{
\tabular{ll}{
Package: \tab ccaPP\cr
Type: \tab Package\cr
Version: \tab 0.2.0\cr
Date: \tab 2012-09-19\cr
Depends: \tab Rcpp (>= 0.9.10), parallel, pcaPP (>= 1.8-1)\cr
Imports: \tab RcppArmadillo (>= 0.2.36)\cr
LinkingTo: \tab Rcpp, RcppArmadillo\cr
Suggests: \tab mvtnorm\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}

Index:
\preformatted{
ccaGrid                 (Robust) CCA via alternating series of grid
                        searches
ccaPP-package           (Robust) canonical correlation analysis via
                        projection pursuit
ccaProj                 (Robust) CCA via projections through the data
                        points
corFunctions            Fast implementations of (robust) correlation
                        estimators
fastMAD                 Fast implementation of the median absolute
                        deviation
fastMedian              Fast implementation of the median
maxCorGrid              (Robust) maximum correlation via alternating
                        series of grid searches
maxCorProj              (Robust) maximum correlation via projections
                        through the data points
permTest                (Robust) permutation test for independence
}
}
\author{
Andreas Alfons [aut, cre],
  David Simcha [ctb]

Maintainer: Andreas Alfons <andreas.alfons@kuleuven.be>
}
\keyword{package}
