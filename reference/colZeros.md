# Number of zeros in columns and rows of large matrices

Produces the same result as applying `sum(x == 0)` to each row or
column.

## Usage

``` r
colZeros(x)

rowZeros(x)
```

## Arguments

- x:

  a [base::matrix](https://rdrr.io/r/base/matrix.html) or
  [Matrix::Matrix](https://rdrr.io/pkg/Matrix/man/Matrix.html) object.

## Examples

``` r
mt <- Matrix::rsparsematrix(100, 100, 0.01)
colZeros(mt)
#>   [1]  97  99 100 100  98  99  99  95  98  98  99  99  99  99  98  98  99  98
#>  [19]  99  99 100 100  99  99 100  99 100 100  98  99  98  99  99  99 100 100
#>  [37] 100  99  97 100 100  99  99 100  98 100 100  99 100  98  99  99 100 100
#>  [55]  99  99 100  98  99 100 100  97  98 100  99 100  96 100  99 100 100 100
#>  [73]  99 100  97  99  97 100 100 100  99  99  98 100 100  99  99  98  97  98
#>  [91] 100 100  98  98  98  99  99  99 100 100
apply(mt, 2, function(x) sum(x == 0)) # the same
#>   [1]  97  99 100 100  98  99  99  95  98  98  99  99  99  99  98  98  99  98
#>  [19]  99  99 100 100  99  99 100  99 100 100  98  99  98  99  99  99 100 100
#>  [37] 100  99  97 100 100  99  99 100  98 100 100  99 100  98  99  99 100 100
#>  [55]  99  99 100  98  99 100 100  97  98 100  99 100  96 100  99 100 100 100
#>  [73]  99 100  97  99  97 100 100 100  99  99  98 100 100  99  99  98  97  98
#>  [91] 100 100  98  98  98  99  99  99 100 100
```
