#' Cross-product of large sparse matrices
#'
#' Compute the (transposed) cross-product of large sparse matrices using the same
#' infrastructure as [simil()] and [dist()].
#' @inheritParams simil
#' @param min_prod the minimum product to be recorded.
#' @export
crossprod <- function(x, y = NULL, min_prod = NULL, sparse = TRUE, digits = 14) {

    if (is.null(min_prod))
        min_prod <- -Inf
    proxy(x, y, margin = 2, method = "product", min_proxy = min_prod,
          sparse = sparse, digits = digits)
}

#' @rdname crossprod
#' @export
tcrossprod <- function(x, y = NULL, min_prod = NULL, sparse = TRUE, digits = 14) {

    if (is.null(min_prod))
        min_prod <- -Inf
    proxy(x, y, margin = 1, method = "product", min_proxy = min_prod,
          sparse = sparse, digits = digits)
}


#' Standard deviation of columns and rows of large matrices
#'
#' Produces the same result as `apply(x, 1, sd)` or `apply(x, 2, sd)`
#' without coercing matrix to dense matrix. Values are not identical to
#' `sd()` because of the floating point precision issue in C++.
#' @param x a [base::matrix] or [Matrix::Matrix] object.
#' @examples
#' mt <- Matrix::rsparsematrix(100, 100, 0.01)
#' colSds(mt)
#' apply(mt, 2, sd) # the same
#' @export
colSds <- function(x) {
    x <- as(as(x, "CsparseMatrix"), "dMatrix")
    structure(cpp_sd(x), names = colnames(x))
}

#' @rdname colSds
#' @export
rowSds <- function(x) {
    x <- as(as(x, "CsparseMatrix"), "dMatrix")
    structure(cpp_sd(t(x)), names = rownames(x))
}

#' Number of zeros in columns and rows of large matrices
#'
#' Produces the same result as applying `sum(x == 0)` to each row or column.
#' @param x a [base::matrix] or [Matrix::Matrix] object.
#' @examples
#' mt <- Matrix::rsparsematrix(100, 100, 0.01)
#' colZeros(mt)
#' apply(mt, 2, function(x) sum(x == 0)) # the same
#' @export
colZeros <- function(x) {
    x <- as(as(x, "CsparseMatrix"), "dMatrix")
    structure(nrow(x) - cpp_nz(x), names = colnames(x))
}

#' @rdname colZeros
#' @export
rowZeros <- function(x) {
    x <- as(as(x, "CsparseMatrix"), "dMatrix")
    structure(ncol(x) - cpp_nz(t(x)), names = rownames(x))
}

#' Maximum/minimum values in columns or rows
#'
#' Produces the same result as `apply(x, 1, max)` or `apply(x, 2, min)`.
#' @param x a [base::matrix] or [Matrix::Matrix] object.
#' @examples
#' mt <- Matrix::rsparsematrix(100, 100, 0.01)
#' colMaxs(mt)
#' apply(mt, 2, max) # the same
#' @export
colMaxs <- function(x) {
    x <- as(as(x, "CsparseMatrix"), "dMatrix")
    structure(cpp_max(x), names = colnames(x))
}

#' @rdname colMaxs
#' @export
rowMaxs <- function(x) {
    x <- as(as(x, "CsparseMatrix"), "dMatrix")
    structure(cpp_max(t(x)), names = rownames(x))
}

#' @rdname colMaxs
#' @export
colMins <- function(x) {
    x <- as(as(x, "CsparseMatrix"), "dMatrix")
    structure(cpp_min(x), names = colnames(x))
}

#' @rdname colMaxs
#' @export
rowMins <- function(x) {
    x <- as(as(x, "CsparseMatrix"), "dMatrix")
    structure(cpp_min(t(x)), names = rownames(x))
}

#' Column or row indices of maximum values
#'
#' Find Column or row indices of maximum values in sparse matrices. Works in a
#' similar way as `base::max.col(ties.method = "first")` but without conversion
#' to dense matrices.
#' @param x a [base::matrix] or [Matrix::Matrix] object.
#' @examples
#' mt <- Matrix::rsparsematrix(100, 100, 0.01)
#' maxCol(mt)
#' max.col(mt, ties.method = "first") # the same
#' @export
maxCol <- function(x) {
    x <- as(as(x, "CsparseMatrix"), "dMatrix")
    structure(cpp_maxrow(t(x)) + 1L, names = rownames(x))
}

#' @rdname maxCol
#' @export
maxRow <- function(x) {
    x <- as(as(x, "CsparseMatrix"), "dMatrix")
    structure(cpp_maxrow(x) + 1L, names = colnames(x))
}

getThreads <- function() {

    # respect other settings
    default <- c("tbb" = as.integer(Sys.getenv("RCPP_PARALLEL_NUM_THREADS")),
                 "omp" = as.integer(Sys.getenv("OMP_THREAD_LIMIT")),
                 "max" = cpp_get_max_thread())
    default <- unname(min(default, na.rm = TRUE))
    suppressWarnings({
        value <- as.integer(getOption("proxyC.threads", default))
    })
    if (length(value) != 1 || is.na(value)) {
        stop("proxyC.threads must be an integer")
    }
    return(value)
}

