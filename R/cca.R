# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' (Robust) CCA via alternating series of grid searches
#' 
#' Perform canoncial correlation analysis via projection pursuit based on 
#' alternating series of grid searches in two-dimensional subspaces of each 
#' data set, with a focus on robust and nonparametric methods.
#' 
#' The algorithm is based on alternating series of grid searches in 
#' two-dimensional subspaces of each data set.  In each grid search, 
#' \code{nGrid} grid points on the unit circle in the corresponding plane are 
#' obtained, and the directions from the center to each of the grid points are 
#' examined.  In the first iteration, equispaced grid points in the interval 
#' \eqn{[-\pi/2, \pi/2)}{[-pi/2, pi/2)} are used.  In each subsequent 
#' iteration, the angles are halved such that the interval 
#' \eqn{[-\pi/4, \pi/4)}{[-pi/4, pi/4)} is used in the second iteration and so 
#' on.  If only one data set is multivariate, the algorithm simplifies 
#' to iterative grid searches in two-dimensional subspaces of the corresponding 
#' data set.
#' 
#' In the basic algorithm, the order of the variables in a series of grid 
#' searches for each of the data sets is determined by the average absolute 
#' correlations with the variables of the respective other data set.  Since 
#' this requires to compute the full \eqn{(p \times q)}{(p x q)} matrix of 
#' absolute correlations, where \eqn{p} denotes the number of variables of 
#' \code{x} and \eqn{q} the number of variables of \code{y}, a faster 
#' modification is available as well.  In this modification, the average 
#' absolute correlations are computed over only a subset of the variables of 
#' the respective other data set.  It is thereby possible to use randomly 
#' selected subsets of variables, or to specify the subsets of variables 
#' directly.
#' 
#' @aliases print.cca
#' 
#' @param x,y  each can be a numeric vector, matrix or data frame.
#' @param k  an integer giving the number of canonical variables to compute.
#' @param method  a character string specifying the correlation functional to 
#' maximize.  Possible values are \code{"spearman"} for the Spearman 
#' correlation, \code{"kendall"} for the Kendall correlation, \code{"quadrant"} 
#' for the quadrant correlation, \code{"M"} for the correlation based on a 
#' bivariate M-estimator of location and scatter with a Huber loss function, or 
#' \code{"pearson"} for the classical Pearson correlation (see 
#' \code{\link{corFunctions}}).
#' @param control  a list of additional arguments to be passed to the specified 
#' correlation functional.  If supplied, this takes precedence over additional 
#' arguments supplied via the \code{\dots} argument.
#' @param nIterations,maxiter  an integer giving the maximum number of 
#' iterations.
#' @param nAlternate,maxalter  an integer giving the maximum number of 
#' alternate series of grid searches in each iteration.
#' @param nGrid,splitcircle  an integer giving the number of equally spaced 
#' grid points on the unit circle to use in each grid search.
#' @param select  optional; either an integer vector of length two or a list 
#' containing two index vectors.  In the first case, the first integer gives 
#' the number of variables of \code{x} to be randomly selected for determining 
#' the order of the variables of \code{y} in the corresponding series of grid 
#' searches, and vice versa for the second integer.  In the latter case, the 
#' first list element gives the indices of the variables of \code{x} to be used 
#' for determining the order of the variables of \code{y}, and vice versa for 
#' the second integer (see \dQuote{Details}).
#' @param tol,zero.tol  a small positive numeric value to be used for 
#' determining convergence.
#' @param seed  optional initial seed for the random number generator (see 
#' \code{\link{.Random.seed}}).  This is only used if \code{select} specifies 
#' the numbers of variables of each data set to be randomly selected for 
#' determining the order of the variables of the respective other data set.
#' @param \dots  additional arguments to be passed to the specified correlation 
#' functional.
#' 
#' @returnClass cca
#' @returnItem cor  a numeric vector giving the canonical correlation 
#' measures.
#' @returnItem A  a numeric matrix in which the columns contain the canonical 
#' vectors for \code{x}.
#' @returnItem B  a numeric matrix in which the columns contain the canonical 
#' vectors for \code{y}.
#' @returnItem call  the matched function call.
#' 
#' @note \code{CCAgrid} is a simple wrapper function for \code{ccaGrid} for 
#' more compatibility with package \pkg{pcaPP} concerning function and argument 
#' names.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{ccaProj}}, \code{\link{corFunctions}}
#' 
#' @examples 
#' ## generate data
#' library("mvtnorm")
#' set.seed(1234)  # for reproducibility
#' p <- 3
#' q <- 2
#' m <- p + q
#' sigma <- 0.5^t(sapply(1:m, function(i, j) abs(i-j), 1:m))
#' xy <- rmvnorm(100, sigma=sigma)
#' x <- xy[, 1:p]
#' y <- xy[, (p+1):m]
#' 
#' ## Spearman correlation
#' ccaGrid(x, y, method = "spearman")
#' ccaGrid(x, y, method = "spearman", consistent = TRUE)
#' 
#' ## Pearson correlation
#' ccaGrid(x, y, method = "pearson")
#' 
#' @keywords multivariate robust
#' 
#' @import Rcpp
#' @import RcppArmadillo
#' @useDynLib ccaPP
#' @export

ccaGrid <- function(x, y, k = 1, 
        method = c("spearman", "kendall", "quadrant", "M", "pearson"), 
        control = list(...), nIterations = 10, nAlternate = 10, nGrid = 25, 
        select = NULL, tol = 1e-06, seed = NULL, ...) {
    ## initializations
    matchedCall <- match.call()
    ## define list of control arguments for algorithm
    nIterations <- as.integer(nIterations)
    nAlternate <- as.integer(nAlternate)
    nGrid <- as.integer(nGrid)
    tol <- as.numeric(tol)
    ppControl <- list(nIterations=nIterations, nAlternate=nAlternate, 
        nGrid=nGrid, select=select, tol=tol)
    ## call workhorse function
    cca <- ccaPP(x, y, k, method=method, corControl=control, 
        algorithm="grid", ppControl=ppControl)
    cca$call <- matchedCall
    cca
}

## wrapper function for more compatibility with package pcaPP
#' @rdname ccaGrid
#' @export

CCAgrid <- function(x, y, k = 1, 
        method = c("spearman", "kendall", "quadrant", "M", "pearson"), 
        maxiter = 10, maxalter = 10, splitcircle = 25, select=NULL, 
        zero.tol = 1e-06, seed = NULL, ...) {
    ## call ccaGrid()
    ccaGrid(x, y, k=k, method=method, nIterations=maxiter, nAlternate=maxalter, 
        nGrid=splitcircle, select=select, tol=zero.tol, seed=seed, ...)
}


#' (Robust) CCA via projections through the data points
#' 
#' Perform canoncial correlation analysis via projection pursuit based on 
#' projections through the data points, with a focus on robust and 
#' nonparametric methods.
#' 
#' First the candidate projection directions are defined for each data set 
#' from the respective center through each data point.  Then the algorithm 
#' scans all \eqn{n^2} possible combinations for the maximum correlation, 
#' where \eqn{n} is the number of observations.
#' 
#' @param x,y  each can be a numeric vector, matrix or data frame.
#' @param k  an integer giving the number of canonical variables to compute.
#' @param method  a character string specifying the correlation functional to 
#' maximize.  Possible values are \code{"spearman"} for the Spearman 
#' correlation, \code{"kendall"} for the Kendall correlation, \code{"quadrant"} 
#' for the quadrant correlation, \code{"M"} for the correlation based on a 
#' bivariate M-estimator of location and scatter with a Huber loss function, or 
#' \code{"pearson"} for the classical Pearson correlation (see 
#' \code{\link{corFunctions}}).
#' @param control  a list of additional arguments to be passed to the specified 
#' correlation functional.  If supplied, this takes precedence over additional 
#' arguments supplied via the \code{\dots} argument.
#' @param useL1Median  a logical indicating whether the \eqn{L_{1}}{L1} medians 
#' should be used as the centers of the data sets (defaults to 
#' \code{TRUE}).  If \code{FALSE}, the columnwise centers are used instead 
#' (columnwise means if \code{method} is \code{"pearson"} and columnwise 
#' medians otherwise).
#' @param \dots  additional arguments to be passed to the specified correlation 
#' functional.
#' 
#' @returnClass cca
#' @returnItem cor  a numeric vector giving the canonical correlation 
#' measures.
#' @returnItem A  a numeric matrix in which the columns contain the canonical 
#' vectors for \code{x}.
#' @returnItem B  a numeric matrix in which the columns contain the canonical 
#' vectors for \code{y}.
#' @returnItem call  the matched function call.
#' 
#' @note \code{CCAproj} is a simple wrapper function for \code{ccaProj} for 
#' more compatibility with package \pkg{pcaPP} concerning function names.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{ccaGrid}} \code{\link{corFunctions}}
#' 
#' @examples 
#' ## generate data
#' library("mvtnorm")
#' set.seed(1234)  # for reproducibility
#' p <- 3
#' q <- 2
#' m <- p + q
#' sigma <- 0.5^t(sapply(1:m, function(i, j) abs(i-j), 1:m))
#' xy <- rmvnorm(100, sigma=sigma)
#' x <- xy[, 1:p]
#' y <- xy[, (p+1):m]
#' 
#' ## Spearman correlation
#' ccaProj(x, y, method = "spearman")
#' ccaProj(x, y, method = "spearman", consistent = TRUE)
#' 
#' ## Pearson correlation
#' ccaProj(x, y, method = "pearson")
#' 
#' @keywords multivariate robust
#' 
#' @import Rcpp
#' @import RcppArmadillo
#' @import pcaPP
#' @useDynLib ccaPP
#' @export

ccaProj <- function(x, y, k = 1, 
        method = c("spearman", "kendall", "quadrant", "M", "pearson"), 
        control = list(...), useL1Median = TRUE, ...) {
    ## initializations
    matchedCall <- match.call()
    ## define list of control arguments for algorithm
    ppControl <- list(useL1Median=isTRUE(useL1Median))
    ## call workhorse function
    cca <- ccaPP(x, y, k, method=method, corControl=control, algorithm="proj", 
        ppControl=ppControl)
    cca$call <- matchedCall
    cca
}

## wrapper function for more compatibility with package pcaPP
#' @rdname ccaProj
#' @export

CCAproj <- function(x, y, k = 1, 
        method = c("spearman", "kendall", "quadrant", "M", "pearson"), 
        useL1Median = TRUE, ...) {
    ## call ccaGrid()
    ccaGrid(x, y, k=k, method=method, useL1Median=useL1Median, ...)
}


## workhorse function
ccaPP <- function(x, y, k = 1, 
        method = c("spearman", "kendall", "quadrant", "M", "pearson"), 
        corControl, algorithm = c("grid", "proj"), ppControl, 
        seed = NULL) {
    ## initializations
    matchedCall <- match.call()
    x <- as.matrix(x)
    y <- as.matrix(y)
    n <- nrow(x)
    if(nrow(y) != n) {
        stop("'x' and 'y' must have the same number of observations")
    }
    p <- ncol(x)
    q <- ncol(y)
    # check number of canonical variables to compute
    k <- rep(as.integer(k), length.out=1)
    if(is.na(k) || k < 0) k <- formals()$k
    k <- min(k, p, q)
    ## prepare the data and call C++ function
    if(n == 0 || p == 0 || q == 0 || k == 0) {
        # zero dimension
        A <- B <- matrix(numeric(), 0, 0)
        cca <- list(cor=NA, A=A, B=B)
    } else {
        # check method and get list of control arguments
        method <- match.arg(method)
        corControl <- getCorControl(method, corControl)
        # additional checks for grid search algorithm
        if(algorithm == "grid") {
            # check subset of variables to be used for determining the order of 
            # the variables from the respective other data set
            select <- ppControl$select
            ppControl$select <- NULL
            if(!is.null(select)) {
                if(is.list(select)) {
                    # make sure select is a list with two index vectors and 
                    # drop invalid indices from each vector
                    select <- rep(select, length.out=2)
                    select <- mapply(function(indices, max) {
                            indices <- as.integer(indices)
                            indices[which(indices > 0 & indices <= max)] - 1
                        }, select, c(p, q))
                    valid <- sapply(select, length) > 0
                    # add the two index vectors to control object
                    if(all(valid)) {
                        ppControl$selectX <- select[[1]]
                        ppControl$selectY <- select[[2]]
                    } else select <- NULL
                } else {
                    # check number of indices to sample
                    select <- rep(as.integer(select), length.out=2)
                    valid <- !is.na(select) & select > 0 & select < c(p, q)
                    if(all(valid)) {
                        # generate index vectors and add them to control object
                        if(!is.null(seed)) set.seed(seed)
                        ppControl$selectX <- sample.int(p, select[1]) - 1
                        ppControl$selectY <- sample.int(q, select[2]) - 1
                    } else select <- NULL
                }
            }
            if(is.null(select)) {
                ppControl$selectX <- ppControl$selectY <- integer()
            }
        }
        # call C++ function
        cca <- .Call("R_ccaPP", R_x=x, R_y=y, R_k=k, R_method=method, 
            R_corControl=corControl, R_algorithm=algorithm, 
            R_ppControl=ppControl, PACKAGE="ccaPP")
        cca$cor <- drop(cca$cor)
    }
    ## assign class and return results
    class(cca) <- "cca"
    cca
}
