#' Create a pattern matrix for masking
#'
#' Create a pattern matrix for `simil()` or `dist()` to enable masking.
#' @param x,y a numeric or character vector to match against each other.
#' @return a sparse logical matrix with `TRUE` for matched pairs.
#' @export
#' @examples
#' mask(c("a", "b", "c", "d", "e"),
#'      c("a", "a", "d", "d", "e", "e"))
#'
mask <- function(x, y) {

    if (!identical(class(x), class(y)))
        stop("x and y must be the same type of vector")
    z <- union(x, y)

    result <- cpp_mask(match(x, z), match(y, z))
    result <- as(result, "lMatrix")
    if (is.character(x) || is.factor(x))
        rownames(result) <- as.character(x)
    if (is.character(y) || is.factor(x))
        colnames(result) <- as.character(y)

    return(result)
}
