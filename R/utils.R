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

#' Number of zeros in columns and rows of large matrices
#'
#' Produce the same result as applying `sum(x == 0)` to each row or column.
#' @param x a [base::matrix] or [Matrix::Matrix] object.
#' @return an integer vector
#' @examples
#' mt <- Matrix::rsparsematrix(100, 100, 0.9)
#' colZeros(mt)
#' apply(mt, 2, function(x) sum(x == 0)) # the same
#' @export
colZeros <- function(x) {
    x <- convertMatrix(x)
    structure(nrow(x) - cpp_nz(x, margin = 2), names = colnames(x))
}

#' @rdname colZeros
#' @export
rowZeros <- function(x) {
    x <- convertMatrix(x)
    structure(ncol(x) - cpp_nz(x, margin = 1), names = rownames(x))
}

#' Standard deviation of columns and rows of large matrices
#'
#' Produce the same result as `apply(x, 2, sd)` and `apply(x, 1, sd)`, respectively.
#' without coercing matrix to dense matrix. Values are not identical to
#' `sd()` because of the floating point precision issue in C++.
#' @param x a [base::matrix] or [Matrix::Matrix] object.
#' @param use_nan if `TRUE`, return `NaN` for columns or rows with all zero.
#'   Only used when `x` is a sparse matrix.
#' @return a numeric vector
#' @examples
#' mt <- Matrix::rsparsematrix(100, 100, 0.9)
#' colSds(mt)
#' apply(mt, 2, sd) # the same
#' @export
colSds <- function(x, use_nan = FALSE) {
    x <- convertMatrix(x)
    use_nan <- check_use_nan(x, use_nan)
    structure(cpp_sd(x, use_nan, margin = 2), names = colnames(x))
}

#' @rdname colSds
#' @export
rowSds <- function(x, use_nan = FALSE) {
    x <- convertMatrix(x)
    use_nan <- check_use_nan(x, use_nan)
    structure(cpp_sd(x, use_nan, margin = 1), names = rownames(x))
}

#' Maximum/minimum values in columns or rows
#'
#' Produce the same result as `apply(x, 2, max)`, `apply(x, 1, max)`,
#' `apply(x, 2, min)` and `apply(x, 1, min)`, respectively, but ignore `NA`.
#' @inheritParams colSds
#' @return a numeric vector
#' @examples
#'   mt <- Matrix::rsparsematrix(100, 100, 0.9)
#'   colMaxs(mt)
#'   apply(mt, 2, max) # the same
#' @export
colMaxs <- function(x, use_nan = FALSE) {
    x <- convertMatrix(x)
    use_nan <- check_use_nan(x, use_nan)
    structure(cpp_max(x, use_nan, margin = 2), names = colnames(x))
}

#' @rdname colMaxs
#' @export
rowMaxs <- function(x, use_nan = FALSE) {
    x <- convertMatrix(x)
    use_nan <- check_use_nan(x, use_nan)
    structure(cpp_max(x, use_nan, margin = 1), names = rownames(x))
}

#' @rdname colMaxs
#' @export
colMins <- function(x, use_nan = FALSE) {
    x <- convertMatrix(x)
    use_nan <- check_use_nan(x, use_nan)
    structure(cpp_min(x, use_nan, margin = 2), names = colnames(x))
}

#' @rdname colMaxs
#' @export
rowMins <- function(x, use_nan = FALSE) {
    x <- convertMatrix(x)
    use_nan <- check_use_nan(x, use_nan)
    structure(cpp_min(x, use_nan, margin = 1), names = rownames(x))
}

#' Column or row indices of maximum values
#'
#' Find column or row indices of maximum values in matrices. Works in a
#' similar way as `base::max.col()` but without conversion of sparse matrices
#' to dense matrices.
#' @param x a [base::matrix] or [Matrix::Matrix] object.
#' @param ties_method specify how to handle ties. If "random", choose an
#'   index randomly within the ties; if "first" or "last" returns the index of the
#'   first or the last tie elements, respectively.
#' @param use_nan if `TRUE`, return `NA` for columns or rows with all zero.
#'   Applies only when `x` is a sparse matrix.
#' @return an integer vector
#' @examples
#' mt <- Matrix::rsparsematrix(100, 100, 0.9)
#' maxCol(mt, ties_method = "first")
#' max.col(mt, ties.method = "first") # the same
#' @export
maxCol <- function(x, ties_method = c("random", "first", "last"),
                   use_nan = FALSE) {
    x <- convertMatrix(x)
    ties_method <- match.arg(ties_method)
    use_nan <- check_use_nan(x, use_nan)

    mode <- match(ties_method, c("random", "first", "last"))
    structure(cpp_maxrow(x, mode, use_nan, margin = 1) + 1L, names = rownames(x))
}

#' @rdname maxCol
#' @export
maxRow <- function(x, ties_method = c("random", "first", "last"),
                   use_nan = FALSE) {
    x <- convertMatrix(x)
    ties_method <- match.arg(ties_method)
    use_nan <- check_use_nan(x, use_nan)

    mode <- match(ties_method, c("random", "first", "last"))
    structure(cpp_maxrow(x, mode, use_nan, margin = 2) + 1L, names = colnames(x))
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

# dense input stays a base double matrix, so it reaches C++ without a copy; everything else
# becomes a dgCMatrix. The boolean methods binarise in C++, together with the transpose for
# margin = 1, so that also takes no extra copy
convertMatrix <- function(x) {
    if (is(x, "denseMatrix") || is.vector(x))
        x <- as.matrix(x)
    if (is.matrix(x)) {
        if (!is.numeric(x) && !is.logical(x))
            stop("x must be a numeric or logical matrix", call. = FALSE)
        if (!is.double(x))
            storage.mode(x) <- "double"
    } else {
        x <- as(as(as(x, "CsparseMatrix"), "generalMatrix"), "dMatrix")
    }
    return(x)
}

check_use_nan <- function(x, use_nan) {
    use_nan <- check_logical(use_nan)
    if (is.matrix(x) && use_nan)
        stop("use_nan applies only x is a sparse matrix", call. = FALSE)
    return(use_nan)
}

check_logical <- function(x) {
    arg <- deparse(substitute(x))
    x <- as.logical(x)
    if (length(x) != 1 || !x %in% c(FALSE, TRUE))
        stop(arg, " must a single logical value", call. = FALSE)
    return(x)
}

check_numeric <- function(x) {
    arg <- deparse(substitute(x))
    x <- as.numeric(x)
    if (length(x) != 1 || is.na(x))
        stop(arg, " must a single numeric value", call. = FALSE)
    return(x)
}

check_integer <- function(x) {
    arg <- deparse(substitute(x))
    x <- as.integer(x)
    if (length(x) != 1 || is.na(x))
        stop(arg, " must a single integer value", call. = FALSE)
    return(x)
}


