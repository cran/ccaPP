# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' @S3method print cca
print.cca <- function(x, ...) {
    # print function call
    if(!is.null(call <- x$call)) {
        cat("\nCall:\n")
        dput(x$call)
    }
    # print maximum correlation
    cat("\nCanonical correlations:\n")
    print(x$cor, ...)
    # return object invisibly
    invisible(x)
}

#' @S3method print maxCor
print.maxCor <- function(x, ...) {
    # print function call
    if(!is.null(call <- x$call)) {
        cat("\nCall:\n")
        dput(x$call)
    }
    # print maximum correlation
    cat("\nMaximum correlation:\n")
    print(x$cor, ...)
    # return object invisibly
    invisible(x)
}

#' @S3method print permTest
print.permTest <- function(x, ...) {
    # print general statement
    cat("\nPermutation test for independence\n\n")
    # print maximum correlation and p-value
    cat(sprintf("r = %f, p-value = %f\n", x$cor0, x$pValue))
    # print number of random permuations
    cat(sprintf("R = %d random permuations\n", x$R))
    # print alternative hypothesis
    cat("Alternative hypothesis: true maximum correlation is not equal to 0\n")
    # return object invisibly
    invisible(x)
}
