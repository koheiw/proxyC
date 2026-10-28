# Maximum/minimum values in columns or rows

Produce the same result as `apply(x, 2, max)`, `apply(x, 1, max)`,
`apply(x, 2, min)` and `apply(x, 1, min)`, respectively, but ignore
`NA`.

## Usage

``` r
colMaxs(x, use_nan = FALSE)

rowMaxs(x, use_nan = FALSE)

colMins(x, use_nan = FALSE)

rowMins(x, use_nan = FALSE)
```

## Arguments

- x:

  a [base::matrix](https://rdrr.io/r/base/matrix.html) or
  [Matrix::Matrix](https://rdrr.io/pkg/Matrix/man/Matrix.html) object.

- use_nan:

  if `TRUE`, return `NaN` for columns or rows with all zero. Only used
  when `x` is a sparse matrix.

## Value

a numeric vector

## Examples

``` r
  mt <- Matrix::rsparsematrix(100, 100, 0.9)
  colMaxs(mt)
#>   [1] 1.7 2.3 2.0 2.3 2.9 2.5 2.4 2.3 1.9 2.3 2.0 1.9 3.4 1.7 2.2 2.4 2.3 2.3
#>  [19] 2.1 3.1 2.2 2.7 1.6 2.0 3.1 2.3 2.4 2.2 2.2 2.1 2.4 2.5 3.5 2.5 2.5 1.8
#>  [37] 2.2 2.9 1.9 2.1 2.8 2.0 2.3 2.6 2.4 2.0 2.3 2.3 2.4 3.1 2.0 2.1 2.7 2.0
#>  [55] 2.3 2.0 2.0 1.8 2.8 2.4 2.5 2.6 2.4 3.6 2.3 2.9 2.8 2.6 2.0 3.1 2.3 2.1
#>  [73] 2.6 2.8 2.2 3.4 3.0 2.3 1.9 1.9 2.1 2.1 2.6 2.1 2.8 2.7 3.2 2.2 2.3 2.4
#>  [91] 2.7 2.1 2.9 2.5 2.5 2.4 2.2 3.2 2.8 2.0
  apply(mt, 2, max) # the same
#>   [1] 1.7 2.3 2.0 2.3 2.9 2.5 2.4 2.3 1.9 2.3 2.0 1.9 3.4 1.7 2.2 2.4 2.3 2.3
#>  [19] 2.1 3.1 2.2 2.7 1.6 2.0 3.1 2.3 2.4 2.2 2.2 2.1 2.4 2.5 3.5 2.5 2.5 1.8
#>  [37] 2.2 2.9 1.9 2.1 2.8 2.0 2.3 2.6 2.4 2.0 2.3 2.3 2.4 3.1 2.0 2.1 2.7 2.0
#>  [55] 2.3 2.0 2.0 1.8 2.8 2.4 2.5 2.6 2.4 3.6 2.3 2.9 2.8 2.6 2.0 3.1 2.3 2.1
#>  [73] 2.6 2.8 2.2 3.4 3.0 2.3 1.9 1.9 2.1 2.1 2.6 2.1 2.8 2.7 3.2 2.2 2.3 2.4
#>  [91] 2.7 2.1 2.9 2.5 2.5 2.4 2.2 3.2 2.8 2.0
```
