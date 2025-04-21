#' Create a pattern matrix for masking
#'
#' Create a pattern matrix for [simil()] or [dist()] to enable masked similarity computation.
#' If the matrix is passed to the function, it computes similarity scores only for cells with `TRUE`.
#' @param x a numeric or character vector matched against each other.
#' @param y a numeric or character vector matched against `x` if provided.
#' @return a sparse logical matrix with `TRUE` for matched pairs.
#' @export
#' @examples
#' mt1 <- Matrix::rsparsematrix(100, 6, 1.0)
#' colnames(mt1) <- c("a", "a", "d", "d", "e", "e")
#' mt2 <- Matrix::rsparsematrix(100, 5, 1.0)
#' colnames(mt2) <- c("a", "b", "c", "d", "e")
#'
#' (msk <- mask(colnames(mt1), colnames(mt2)))
#' simil(mt1, mt2, margin = 2, mask = msk, drop0 = TRUE)
mask <- function(x, y = NULL) {

    if (is.null(y))
        y <- x
    if (!identical(class(x), class(y)))
        stop("x and y must be the same type of vectors")
    z <- union(x, y)

    result <- cpp_mask(match(x, z), match(y, z), getThreads())
    result <- as(result, "lMatrix")
    if (is.character(x) || is.factor(x))
        rownames(result) <- as.character(x)
    if (is.character(y) || is.factor(x))
        colnames(result) <- as.character(y)

    return(result)
}
