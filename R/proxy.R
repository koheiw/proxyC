#' Compute similarity/distance between rows or columns of large matrices
#'
#' Fast similarity/distance computation function for large sparse matrices. You
#' can floor small similarity value to to save computation time and storage
#' space by an arbitrary threshold (\code{min_simil}) or rank (\code{rank}). You
#' can specify the number of threads for parallel computing via
#' `options(proxyC.threads)`.
#'
#' @param x \link{matrix} or \link{Matrix} object. Dense matrices are covered to
#'   the \link{CsparseMatrix-class} internally.
#' @param y if a \link{matrix} or \link{Matrix} object is provided, proximity
#'   between documents or features in \code{x} and \code{y} is computed.
#' @param margin integer indicating margin of similarity/distance computation. 1
#'   indicates rows or 2 indicates columns.
#' @param method method to compute similarity or distance
#' @param min_simil the minimum similarity value to be recorded.
#' @param rank an integer value specifying top-n most similarity values to be
#'   recorded.
#' @param p weight for Minkowski distance
#' @param drop0 if \code{TRUE}, removes zero values to make the
#'   similarity/distance matrix sparse. It has no effect when `dense = TRUE`.
#' @param diag if \code{TRUE}, only compute diagonal elements of the
#'   similarity/distance matrix; useful when comparing corresponding rows or
#'   columns of `x` and `y`.
#' @param use_nan if `TRUE`, returns `NaN` if the standard deviation of a vector
#'   is zero when `method` is "correlation"; if all the values are zero in a
#'   vector when `method` is "cosine", "chisquared", "kullback", "jeffreys" or
#'   "jensen". Note that use of `NaN` makes the similarity/distance matrix
#'   denser and therefore larger in RAM. If `FALSE`, return zero in same use
#'   situations as above. If `NULL`, will also return zero but also generate a
#'   warning (default).
#' @param sparse if `TRUE`, returns \link{sparseMatrix} object. When neither
#'   `min_simil` nor `rank` is used, dense matrices require less space in RAM.
#' @param digits determines rounding of small values towards zero. Use primarily
#'   to correct floating point errors. Rounding is performed in C++ in a similar
#'   way as \link{zapsmall}.
#' @details ## Available Methods
#'
#'   Similarity:
#' \itemize{
#'   \item `cosine`: cosine similarity
#'   \item `correlation`: Pearson's correlation
#'   \item `jaccard`: Jaccard coefficient
#'   \item `ejaccard`: the real value version of `jaccard`
#'   \item `fjaccard`: Fuzzy Jaccard coefficient
#'   \item `dice`: Dice coefficient
#'   \item `edice`: the real value version of `dice`
#'   \item `hamann`: Hamann similarity
#'   \item `faith`: Faith similarity
#'   \item `simple matching`: the percentage of common elements
#' }
#'   Distance:
#' \itemize{
#'   \item `euclidean`: Euclidean distance
#'   \item `chisquared`: chi-squared distance
#'   \item `kullback`: Kullback–Leibler divergence
#'   \item `jeffreys`: Jeffreys divergence
#'   \item `jensen`: Jensen–Shannon divergence
#'   \item `manhattan`: Manhattan distance
#'   \item `maximum`: the largest difference between values
#'   \item `canberra`: Canberra distance
#'   \item `minkowski`: Minkowski distance
#'   \item `hamming`: Hamming distance
#' }
#'   See the vignette for how the similarity and distance are computed:
#'   `vignette("measures", package = "proxyC")`
#'
#'   ## Parallel Computing
#'
#'   It performs parallel computing using Intel oneAPI Threads Building Blocks.
#'   The number of threads for parallel computing should be specified via
#'   `options(proxyC.threads)` before calling the functions. If the value is -1,
#'   all the available threads will be used. Unless the option is used, the
#'   number of threads will be limited by the environmental variables
#'   (`OMP_THREAD_LIMIT` or `RCPP_PARALLEL_NUM_THREADS`) to comply with CRAN
#'   policy and offer backward compatibility.
#'
#' @import methods Matrix
#' @seealso zapsmall
#' @export
#' @examples
#' mt <- Matrix::rsparsematrix(100, 100, 0.01)
#' simil(mt, method = "cosine")[1:5, 1:5]
simil <- function(x, y = NULL, margin = 1,
                  method = c("cosine", "correlation", "jaccard", "ejaccard", "fjaccard",
                             "dice", "edice", "hamann", "faith", "simple matching"),
                  min_simil = NULL, rank = NULL, drop0 = FALSE, diag = FALSE,
                  use_nan = NULL, sparse = TRUE, digits = 14) {

    method[method == "hamman"] <- "hamann" # for transition
    method <- match.arg(method)
    proxy(x, y, margin, method, min_proxy = min_simil, rank = rank, drop0 = drop0,
          diag = diag, use_nan = use_nan, sparse = sparse, digits = digits)

}

#' @rdname simil
#' @param smooth adds a  fixed value to all the cells to avoid division by zero.
#'   Only used when `method` is "chisquared", "kullback", "jeffreys" or "jensen".

#' @export
#' @examples
#' mt <- Matrix::rsparsematrix(100, 100, 0.01)
#' dist(mt, method = "euclidean")[1:5, 1:5]
dist <- function(x, y = NULL, margin = 1,
                 method = c("euclidean", "chisquared", "kullback", "jeffreys", "jensen",
                            "manhattan", "maximum", "canberra", "minkowski", "hamming"),
                 p = 2, smooth = 0, drop0 = FALSE, diag = FALSE, use_nan = NULL,
                 sparse = TRUE, digits = 14) {

    method <- match.arg(method)
    proxy(x, y, margin, method, p = p, smooth = smooth, drop0 = drop0,
          diag = diag, use_nan = use_nan, sparse = sparse, digits = digits)
}

#' @import Rcpp
#' @useDynLib proxyC
proxy <- function(x, y = NULL, margin = 1,
                  method = c("cosine", "correlation", "jaccard", "ejaccard", "fjaccard",
                             "dice", "edice", "hamann", "simple matching", "faith",
                             "euclidean", "chisquared", "kullback", "jeffreys", "jensen",
                             "manhattan", "maximum", "canberra", "minkowski", "hamming"),
                  p = 2, smooth = 0, min_proxy = NULL, rank = NULL, drop0 = FALSE,
                  diag = FALSE, use_nan = NULL, sparse = TRUE, digits = 14) {

    method[method == "hamman"] <- "hamann" # for transition
    method <- match.arg(method)
    #x <- as(as(x, "CsparseMatrix"), "dgCMatrix")
    x <- as(as(as(x, "CsparseMatrix"), "generalMatrix"), "dMatrix") # for Matrix v1.4-2 or later
    symm <- is.null(y)

    if (is.null(y)) {
        y <- x
    } else {
        #y <- as(as(y, "CsparseMatrix"), "dgCMatrix")
        y <- as(as(as(y, "CsparseMatrix"), "generalMatrix"), "dMatrix") # for Matrix v1.4-2 or later
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
    if (is.null(use_nan)) {
        if (method == "correlation" && (any(colSds(x) == 0) || any(colSds(y) == 0))) {
            warning(paste0(
                "x or y has vectors with zero standard deviation; ",
                "consider setting use_nan = TRUE to set these values to NaN ",
                "or use_nan = FALSE to suppress this warning"), call. = FALSE)
        } else if (
            method %in% c("cosine", "kullback", "chisquared", "jeffreys", "jensen") &&
            (any(colZeros(x) == nrow(x)) || any(colZeros(y) == nrow(y)))) {
            warning(paste0(
                "x or y has vectors with all zero; ",
                "consider setting use_nan = TRUE to set these values to NaN ",
                "or use_nan = FALSE to suppress this warning"), call. = FALSE)
        }
        use_nan <- FALSE
    }
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
    } else if (method == "hamann") {
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
        x <- as(as(x, "lMatrix"), "dMatrix")
        y <- as(as(y, "lMatrix"), "dMatrix")
    }
    if (method %in% c("cosine", "correlation", "euclidean") && !diag) {
        result <- cpp_linear(
            mt1 = x,
            mt2 = y,
            method = match(method, c("cosine", "correlation", "euclidean")),
            rank = rank,
            limit = min_proxy,
            symm = symm,
            drop0 = drop0,
            use_nan = use_nan,
            sparse = sparse,
            digits = digits,
            thread = getThreads()
        )
    } else {
        result <- cpp_pair(
            mt1 = x,
            mt2 = y,
            method = match(method, c("cosine", "correlation", "ejaccard", "edice",
                                     "hamann", "simple matching", "faith",
                                     "euclidean", "chisquared", "kullback", "manhattan",
                                     "maximum", "canberra", "minkowski", "hamming",
                                     "jeffreys", "jensen", "fjaccard")),
            rank = rank,
            limit = min_proxy,
            weight = weight,
            smooth = smooth,
            symm = symm,
            diag = diag,
            drop0 = drop0,
            use_nan = use_nan,
            sparse = sparse,
            digits = digits,
            thread = getThreads()
        )
    }
    if (diag)
        result <- as(as(result, "diagonalMatrix"), "ddiMatrix")
    dimnames(result) <- list(colnames(x), colnames(y))
    return(result)
}

#' Standard deviation of columns and rows of large matrices
#'
#' Produces the same result as \code{apply(x, 1, sd)} or \code{apply(x, 2, sd)}
#' without coercing matrix to dense matrix. Values are not identical to
#' \code{sd} because of the floating point precision issue in C++.
#' @param x \link{matrix} or \link{Matrix} object
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

#' Count number of zeros in columns and rows of large matrices
#'
#' Produces the same result as applying \code{sum(x == 0)} to each row or column.
#' @param x \link{matrix} or \link{Matrix} object
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

