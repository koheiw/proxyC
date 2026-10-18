# Create a pattern matrix for masking

Create a pattern matrix for
[`simil()`](https://koheiw.github.io/proxyC/reference/simil.md) or
[`dist()`](https://koheiw.github.io/proxyC/reference/simil.md) to enable
masked similarity computation. If the matrix is passed to the function,
it computes similarity scores only for cells with `TRUE`.

## Usage

``` r
mask(x, y = NULL)
```

## Arguments

- x:

  a numeric or character vector matched against each other.

- y:

  a numeric or character vector matched against `x` if provided.

## Value

a sparse logical matrix with `TRUE` for matched pairs.

## Examples

``` r
mt1 <- Matrix::rsparsematrix(100, 6, 1.0)
colnames(mt1) <- c("a", "a", "d", "d", "e", "e")
mt2 <- Matrix::rsparsematrix(100, 5, 1.0)
colnames(mt2) <- c("a", "b", "c", "d", "e")

(msk <- mask(colnames(mt1), colnames(mt2)))
#> 6 x 5 sparse Matrix of class "lgTMatrix"
#>   a b c d e
#> a | . . . .
#> a | . . . .
#> d . . . | .
#> d . . . | .
#> e . . . . |
#> e . . . . |
simil(mt1, mt2, margin = 2, mask = msk, drop0 = TRUE)
#> 6 x 5 sparse Matrix of class "dgTMatrix"
#>             a b c           d           e
#> a  0.14040326 . . .            .         
#> a -0.02487187 . . .            .         
#> d  .          . . 0.014260459  .         
#> d  .          . . 0.007621309  .         
#> e  .          . . .           -0.06406566
#> e  .          . . .           -0.08833690
```
