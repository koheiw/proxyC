#' Compute the cross product of large sparse matrices
#'
#' Compute the (transposed) cross product of large sparse matrices using the same
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
