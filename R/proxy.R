#' Compute similiarty/distance between raws or columns of large matrices

#' @param x a \link{matrix} or \link{Matrix} object
#' @param y if a \link{matrix} or \link{Matrix} object is provided, proximity
#'   between documents or features in \code{x} and \code{y} is computed.
#' @param method to compute similarity or distance
#' @param min_proxy the minimum proximity value to be recoded.
#' @param rank an integer value specifying top-n most proximity values to be
#'   recorded.
#' @param p weight for minkowski distance
#' @import methods Matrix
#' @export
#' @examples
#' mt <- Matrix::rsparsematrix(100, 100, 0.01)
#' simil(mt, method = "cosine")
simil <- function(x, y = NULL, margin = 1,
                  method = c("cosine", "correlation", "jaccard", "ejaccard",
                             "dice", "edice", "hamman", "simple matching", "faith"),
                  min_proxy = NULL, rank = NULL) {

    UseMethod("simil")

}

#' @rdname simil
#' @export
#' @examples
#' dist(mt, method = "euclidean")
dist <- function(x, y = NULL, margin = 1,
                 method = c("euclidean", "chisquared", "hamming", "kullback",
                            "manhattan", "maximum", "canberra", "minkowski"),
                 p = 2) {

    method <- match.arg(method)
    UseMethod("dist")

}

#' @export
simil.matrix <- function(x, y = NULL, margin = 1,
                         method = c("cosine", "correlation", "jaccard", "ejaccard",
                                    "dice", "edice", "hamman", "simple matching", "faith"),
                         min_proxy = NULL, rank = NULL) {
    simil(as(x, "dgCMatrix"), as(y, "dgCMatrix"), min_proxy, rank)
}

#' @export
dist.matrix <- function(x, y = NULL, margin = 1,
                        method = c("euclidean", "chisquared", "hamming", "kullback",
                                   "manhattan", "maximum", "canberra", "minkowski"),
                        p = 2) {
    dist(as(x, "dgCMatrix"), as(y, "dgCMatrix"), p)
}

#' @export
#' @import Matrix
simil.Matrix <- function(x, y = NULL, margin = 1,
                         method = c("cosine", "correlation", "jaccard", "ejaccard",
                                    "dice", "edice", "hamman", "simple matching", "faith"),
                         min_proxy = NULL, rank = NULL) {
    method <- match.arg(method)
    proxy(x, y, margin, method, min_proxy = min_proxy, rank = rank)
}

#' @export
#' @import Matrix
dist.Matrix <- function(x, y = NULL, margin = 1,
                        method = c("euclidean", "chisquared", "hamming", "kullback",
                                   "manhattan", "maximum", "canberra", "minkowski"), p = 2) {
    method <- match.arg(method)
    proxy(x, y, margin, method, p = p)
}

#' @import RcppParallel Matrix
proxy <- function(x, y = NULL, margin = 1,
                  method = c("cosine", "correlation", "jaccard", "ejaccard",
                             "dice", "edice", "hamman", "simple matching", "faith",
                             "euclidean", "chisquared", "hamming", "kullback",
                             "manhattan", "maximum", "canberra", "minkowski"),
                  p = 2, min_proxy = NULL, rank = NULL) {

    method <- match.arg(method)

    if (is.null(y))
        y <- x
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
        x <- as(x, "lgCMatrix")
        y <- as(y, "lgCMatrix")
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
