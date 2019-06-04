#' Compute similarity/distance between rows or columns of large matrices
#'
#' Fast similarity/distance computation function for large sparse matrices. You
#' can floor small similarity value to to save computation time and storage
#' space by an arbitrary threashold (\code{min_simil}) or rank (\code{rank}).
#' Please increase the numbner of threads for better perfromance using
#' \code{\link[RcppParallel]{setThreadOptions}}.
#'
#' @param x \link{Matrix} object
#' @param y if a \link{matrix} or \link{Matrix} object is provided, proximity
#'   between documents or features in \code{x} and \code{y} is computed.
#' @param margin integer indicating margin of similarity/distance computation. 1
#'   indicates rows or 2 indicates columns.
#' @param method method to compute similarity or distance
#' @param min_simil the minimum similarity value to be recoded.
#' @param rank an integer value specifying top-n most similarity values to be
#'   recorded.
#' @param p weight for minkowski distance
#' @param digits determines rounding of small values towards zero. Use
#'   primarily to correct rounding errors in C++. See \link{zapsmall}.
#' @import methods Matrix
#' @importFrom RcppParallel RcppParallelLibs
#' @seealso zapsmall
#' @export
#' @examples
#' mt <- Matrix::rsparsematrix(100, 100, 0.01)
#' simil(mt, method = "cosine")[1:5, 1:5]
simil <- function(x, y = NULL, margin = 1,
                  method = c("cosine", "correlation", "jaccard", "ejaccard",
                             "dice", "edice", "hamman", "simple matching", "faith"),
                  min_simil = NULL, rank = NULL, digits = 14) {

    method <- match.arg(method)
    proxy(x, y, margin, method, min_proxy = min_simil, rank = rank, digits = digits)

}

#' @rdname simil
#' @export
#' @examples
#' mt <- Matrix::rsparsematrix(100, 100, 0.01)
#' dist(mt, method = "euclidean")[1:5, 1:5]
dist <- function(x, y = NULL, margin = 1,
                 method = c("euclidean", "chisquared", "hamming", "kullback",
                            "manhattan", "maximum", "canberra", "minkowski"),
                 p = 2, digits = 14) {

    method <- match.arg(method)
    proxy(x, y, margin, method, p = p)

}

#' @import Rcpp
#' @useDynLib proxyC
proxy <- function(x, y = NULL, margin = 1,
                  method = c("cosine", "correlation", "jaccard", "ejaccard",
                             "dice", "edice", "hamman", "simple matching", "faith",
                             "euclidean", "chisquared", "hamming", "kullback",
                             "manhattan", "maximum", "canberra", "minkowski"),
                  p = 2, min_proxy = NULL, rank = NULL, digits = 14) {

    method <- match.arg(method)
    if(is(x, 'sparseMatrix')) {
        x <- as(x, "dgCMatrix")
    } else {
        stop("x must be a sparseMatrix")
    }

    symm <- is.null(y)

    if (is.null(y)) {
        y <- x
    } else {
        if(is(y, 'sparseMatrix')) {
            y <- as(y, "dgCMatrix")
        } else {
            stop("y must be a sparseMatrix")
        }
    }
    if (!margin %in% c(1, 2))
        stop("Matrgin must be 1 (row) or 2 (column)")
    if (margin == 1) {
        if (ncol(x) != ncol(y))
            stop("x and y must have the same number of columns")
        x <- t(x)
        y <- t(y)
    } else {
        if (nrow(x) != nrow(y))
            stop("x and y must have the same number of rows")
    }
    if (is.null(min_proxy))
        min_proxy <- -1.0
    if (is.null(rank))
        rank <- ncol(x)
    if (rank < 1)
        stop("rank must be great than or equal to 1")

    boolean <- FALSE
    weight <- 1
    if (method == "jaccard") {
        boolean <- TRUE
        method <- "ejaccard"
    } else if (method == "ejaccard") {
        weight <- 2
    } else if (method == "dice") {
        boolean <- TRUE
        method <- "edice"
    } else if (method == "edice") {
        weight <- 2
    } else if (method == "hamman") {
        boolean <- TRUE
    } else if (method == "faith") {
        boolean <- TRUE
    } else if (method == "simple matching") {
        boolean <- TRUE
    } else if (method == "minkowski") {
        if (p <= 0)
            stop("p must be greater than zero")
        weight <- p
    }
    if (boolean) {
        x <- as(as(x, "lgCMatrix"), "dgCMatrix")
        y <- as(as(y, "lgCMatrix"), "dgCMatrix")
    }
    if (method %in% c("cosine", "correlation", "euclidean")) {
        result <- cpp_linear(
            mt1 = x,
            mt2 = y,
            method = match(method, c("cosine", "correlation", "euclidean")),
            rank = rank,
            limit = min_proxy,
            symm = symm
        )
    } else {
        result <- cpp_pair(
            mt1 = x,
            mt2 = y,
            method = match(method, c("ejaccard", "edice", "hamman", "simple matching",
                                     "faith", "chisquared", "kullback", "manhattan",
                                     "maximum", "canberra", "minkowski")),
            rank = rank,
            limit = min_proxy,
            weight = weight,
            symm = symm
        )
    }
    result@x <- zapsmall(result@x, digits)
    dimnames(result) <- list(colnames(x), colnames(y))
    return(result)
}

#' Standard deviasion of columns and rows in sparse matrix
#'
#' Produces the same result as \code{apply(x, 1, sd)} or \code{apply(x, 2, sd)}
#' as without coecirng matrix to dense matrix. Values are not identical to
#' \code{sd} because of the floating point precision issue in C++.
#' @param x \link{Matrix} object
#' @examples
#' mt <- Matrix::rsparsematrix(100, 100, 0.01)
#' colSds(mt)
#' apply(mt, 2, sd) # the same
#' @export
colSds <- function(x) {
    if(is(x, 'sparseMatrix')) {
        x <- as(x, "dgCMatrix")
    } else {
        stop("x must be a sparseMatrix")
    }
    cpp_sd(x)
}

#' @rdname colSds
#' @export
rowSds <- function(x) {
    if(is(x, 'sparseMatrix')) {
        x <- as(x, "dgCMatrix")
    } else {
        stop("x must be a sparseMatrix")
    }
    cpp_sd(t(x))
}

#' Count number of zeros in columns and rows in sparse matrix
#'
#' Produces the same result as applying \code{sum(x == 0)} to each row or column.
#' @param x \link{Matrix} object
#' @examples
#' mt <- Matrix::rsparsematrix(100, 100, 0.01)
#' colZeros(mt)
#' apply(mt, 2, function(x) sum(x == 0)) # the same
#' @export
colZeros <- function(x) {
    if(is(x, 'sparseMatrix')) {
        x <- as(x, "dgCMatrix")
    } else {
        stop("x must be a sparseMatrix")
    }
    nrow(x) - cpp_nz(x)
}

#' @rdname colZeros
#' @export
rowZeros <- function(x) {
    if(is(x, 'sparseMatrix')) {
        x <- as(x, "dgCMatrix")
    } else {
        stop("x must be a sparseMatrix")
    }
    ncol(x) - cpp_nz(t(x))
}
