#' Cross-product of large sparse matrices
#'
#' Compute the (transposed) cross-product of large sparse matrices using the same
#' infrastructure as [simil()] and [dist()].
#' @inheritParams simil
#' @param min_prod the minimum product to be recorded.
#' @export
crossprod <- function(x, y = NULL, min_prod = NULL, digits = 14) {

    if (is.null(min_prod))
        min_prod <- -Inf
    proxy(x, y, margin = 2, method = "product", min_proxy = min_prod,
          digits = digits)
}

#' @rdname crossprod
#' @export
tcrossprod <- function(x, y = NULL, min_prod = NULL, digits = 14) {

    if (is.null(min_prod))
        min_prod <- -Inf
    proxy(x, y, margin = 1, method = "product", min_proxy = min_prod,
          digits = digits)
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
    result <- cpp_sd(x)
    names(result) <- colnames(x)
    return(result)
}

#' @rdname colSds
#' @export
rowSds <- function(x) {
    x <- as(as(x, "CsparseMatrix"), "dMatrix")
    result <- cpp_sd(t(x))
    names(result) <- rownames(x)
    return(result)
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
    result <- nrow(x) - cpp_nz(x)
    names(result) <- colnames(x)
    return(result)
}

#' @rdname colZeros
#' @export
rowZeros <- function(x) {
    x <- as(as(x, "CsparseMatrix"), "dMatrix")
    result <- ncol(x) - cpp_nz(t(x))
    names(result) <- rownames(x)
    return(result)
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

