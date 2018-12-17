# Compute similiarty/distance between raws or columns of large matrices

#' @export
#' @examples
#' mt <- Matrix::rsparsematrix(100, 100, 0.01)
#' simil(mt)
simil <- function(x, y = NULL, margin = 1,
                  method = c("cosine", "correlation", "jaccard", "ejaccard",
                             "dice", "edice", "hamman", "simple matching", "faith"),
                  min_proxy = NULL, rank = NULL) {

    UseMethod("simil")

}

#' @export
dist <- function(x, y = NULL, margin = 1,
                 method = c("euclidean", "chisquared", "hamming", "kullback",
                            "manhattan", "maximum", "canberra", "minkowski"),
                 p = 2, min_proxy = NULL, rank = NULL) {

    method <- match.arg(method)
    UseMethod("simil")

}

#' @export
#' @import Matrix
simil.matrix <- function(x, y = NULL, margin = 1, ...) {
    simil(as(x, "dgCMatrix"), as(y, "dgCMatrix"), ...)
}

#' @export
#' @import Matrix
dist.matrix <- function(x, y = NULL, margin = 1, ...) {
    simil(as(x, "dgCMatrix"), as(y, "dgCMatrix"), ...)
}

#' @export
#' @import Matrix
simil.Matrix <- function(x, y = NULL, margin = 1,
                         method = c("cosine", "correlation", "jaccard", "ejaccard",
                                    "dice", "edice", "hamman", "simple matching", "faith"), ...) {
    method <- match.arg(method)
    proxy(x, y, margin, method, ...)
}

#' @export
#' @import Matrix
dist.Matrix <- function(x, y = NULL, margin = 1,
                        method = c("euclidean", "chisquared", "hamming", "kullback",
                                   "manhattan", "maximum", "canberra", "minkowski"), ...) {
    method <- match.arg(method)
    proxy(x, y, margin, method, ...)
}

#' [Experimental] Compute document/feature proximity
#'
#' This is an underlying function for \code{textstat_dist} and
#' \code{textstat_simil} but returns \code{TsparseMatrix}.
#' @keywords internal
#' @param y if a \link{dfm} object is provided, proximity between documents or
#'   features in \code{x} and \code{y} is computed.
#' @inheritParams textstat_dist
#' @param min_proxy the minimum proximity value to be recoded.
#' @param rank an integer value specifying top-n most proximity values to be
#'   recorded.
proxy <- function(x, y = NULL, margin = 1,
                  method = c("cosine", "correlation", "jaccard", "ejaccard",
                             "dice", "edice", "hamman", "simple matching", "faith",
                             "euclidean", "chisquared", "hamming", "kullback",
                             "manhattan", "maximum", "canberra", "minkowski"),
                  p = 2, min_proxy = NULL, rank = NULL) {

    if (is.null(y)) {
        y <- x
    } else {
        if (!is.dfm(y))
            stop("y must be a dfm")
        y <- as.dfm(y)
        if (!sum(y)) stop(message_error("dfm_empty"))
    }

    if (margin %in% c(1, 2))
        stop("Matrgin must be 1 (row) or 2 (column)")
    method <- match.arg(method)

    if (margin == 1) {
        f <- union(colnames(x), colnames(y))
        x <- t(pad_dfm(x, f))
        y <- t(pad_dfm(y, f))
    } else {
        if (!identical(rownames(x), rownames(y)))
            stop("x and y must contain the same documents")
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
        x <- dfm_weight(x, "boolean")
        y <- dfm_weight(y, "boolean")
    }
    if (method %in% c("cosine", "correlation", "euclidean")) {
        result <- prxc_linear(x, y,
                              match(method, c("cosine", "correlation", "euclidean")),
                              rank, min_proxy)
    } else {
        result <- prxc_pair(x, y,
                            match(method, c("ejaccard", "edice", "hamman", "simple matching",
                                            "faith", "chisquared", "kullback", "manhattan",
                                            "maximum", "canberra", "minkowski")),
                            rank, min_proxy, weight)
    }

    dimnames(result) <- list(colnames(x), colnames(y))
    return(result)
    # return(as(result, "CsparseMatrix"))
}
