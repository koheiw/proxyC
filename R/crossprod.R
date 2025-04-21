#' Compute the cross product of sparse matrices
#' @inheritParams simil
#' @param min_prod the minimum product to be recorded.
#' @export
crossprod <- function(x, y = NULL, min_prod = NULL, digits = 14) {

    if (is.null(min_prod))
        min_prod <- -Inf
    proxy(x, y, margin = 2, method = "product", min_proxy = min_prod,
          digits = digits)
}

