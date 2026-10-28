# Number of zeros in columns and rows of large matrices

Produce the same result as applying `sum(x == 0)` to each row or column.

## Usage

``` r
colZeros(x)

rowZeros(x)
```

## Arguments

- x:

  a [base::matrix](https://rdrr.io/r/base/matrix.html) or
  [Matrix::Matrix](https://rdrr.io/pkg/Matrix/man/Matrix.html) object.

## Value

an integer vector

## Examples

``` r
mt <- Matrix::rsparsematrix(100, 100, 0.9)
colZeros(mt)
#>   [1] 10 15 12 12 10  9  8 14 10 12 10  7  8 11 11 14  7 13  8 12 10 12 14 11  7
#>  [26]  7  5 11  9 13 11 14 10 11  7 10 11  9  6 16  7  9  8 11  6  6  9  8 16 12
#>  [51] 10 10  7 14 12 17  9 10 14 12 10  7 11 10 12 12  6 13 11  9  6 10 10 11 11
#>  [76]  5  7 10 12  7 11  8  6  9  8  9  7 12 10 11  4  9 13 10  8 10  8 10 12 11
apply(mt, 2, function(x) sum(x == 0)) # the same
#>   [1] 10 15 12 12 10  9  8 14 10 12 10  7  8 11 11 14  7 13  8 12 10 12 14 11  7
#>  [26]  7  5 11  9 13 11 14 10 11  7 10 11  9  6 16  7  9  8 11  6  6  9  8 16 12
#>  [51] 10 10  7 14 12 17  9 10 14 12 10  7 11 10 12 12  6 13 11  9  6 10 10 11 11
#>  [76]  5  7 10 12  7 11  8  6  9  8  9  7 12 10 11  4  9 13 10  8 10  8 10 12 11
```
