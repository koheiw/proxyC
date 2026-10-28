# Column or row indices of maximum values

Find column or row indices of maximum values in matrices. Works in a
similar way as [`base::max.col()`](https://rdrr.io/r/base/maxCol.html)
but without conversion of sparse matrices to dense matrices.

## Usage

``` r
maxCol(x, ties_method = c("random", "first", "last"), use_nan = FALSE)

maxRow(x, ties_method = c("random", "first", "last"), use_nan = FALSE)
```

## Arguments

- x:

  a [base::matrix](https://rdrr.io/r/base/matrix.html) or
  [Matrix::Matrix](https://rdrr.io/pkg/Matrix/man/Matrix.html) object.

- ties_method:

  specify how to handle ties. If "random", choose an index randomly
  within the ties; if "first" or "last" returns the index of the first
  or the last tie elements, respectively.

- use_nan:

  if `TRUE`, return `NA` for columns or rows with all zero. Applies only
  when `x` is a sparse matrix.

## Value

an integer vector

## Examples

``` r
mt <- Matrix::rsparsematrix(100, 100, 0.9)
maxCol(mt, ties_method = "first")
#>   [1] 59 17 25 34 25 64  5 76 30 97 10 44 18 98 54 99 21 47 32 58 58 23  8 73 25
#>  [26] 28 95  1 22 79 32  4 15 70  6 71 11 37 20 65 52 57 47  9 73 40 49 20 54 21
#>  [51]  8 82 31 54 75 91 67 53 33 77 90 29 85 97 97 93 96 55 65 32 15 78 12 92  8
#>  [76] 30 59 53 94 13 81 59 89 93 96 93 69 33 95 14 38 98 37  6 15 88 18 43 17 77
max.col(mt, ties.method = "first") # the same
#>   [1] 59 17 25 34 25 64  5 76 30 97 10 44 18 98 54 99 21 47 32 58 58 23  8 73 25
#>  [26] 28 95  1 22 79 32  4 15 70  6 71 11 37 20 65 52 57 47  9 73 40 49 20 54 21
#>  [51]  8 82 31 54 75 91 67 53 33 77 90 29 85 97 97 93 96 55 65 32 15 78 12 92  8
#>  [76] 30 59 53 94 13 81 59 89 93 96 93 69 33 95 14 38 98 37  6 15 88 18 43 17 77
```
