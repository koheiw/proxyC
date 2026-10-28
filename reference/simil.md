# Compute similarity/distance between rows or columns of large matrices

Fast similarity/distance computation function for large matrices using
multiple processors. To save computation time and storage, it can also
apply arbitrary threshold (`min_simil`) or rank (`rank`) to the
similarity/distance scores and return them in a sparse matrix.

## Usage

``` r
simil(
  x,
  y = NULL,
  margin = 1,
  method = c("cosine", "correlation", "dice", "edice", "jaccard", "ejaccard", "fjaccard",
    "hamann", "faith", "simple matching"),
  mask = NULL,
  min_simil = NULL,
  rank = NULL,
  drop0 = FALSE,
  diag = FALSE,
  use_nan = NULL,
  sparse = TRUE,
  digits = 14
)

dist(
  x,
  y = NULL,
  margin = 1,
  method = c("euclidean", "chisquared", "kullback", "jeffreys", "jensen", "manhattan",
    "maximum", "canberra", "minkowski", "hamming"),
  mask = NULL,
  p = 2,
  smooth = 0,
  drop0 = FALSE,
  diag = FALSE,
  use_nan = NULL,
  sparse = TRUE,
  digits = 14
)
```

## Arguments

- x:

  a [base::matrix](https://rdrr.io/r/base/matrix.html) or
  [Matrix::Matrix](https://rdrr.io/pkg/Matrix/man/Matrix.html) object.
  Only sparse matrices are covered to the
  [Matrix::CsparseMatrix](https://rdrr.io/pkg/Matrix/man/CsparseMatrix-class.html)
  internally.

- y:

  if a [base::matrix](https://rdrr.io/r/base/matrix.html) or
  [Matrix::Matrix](https://rdrr.io/pkg/Matrix/man/Matrix.html) object is
  provided, proximity between documents or features in `x` and `y` is
  computed.

- margin:

  integer indicating margin of similarity/distance computation. 1
  indicates rows or 2 indicates columns.

- method:

  method to compute similarity or distance

- mask:

  a pattern matrix created using
  [`mask()`](https://koheiw.github.io/proxyC/reference/mask.md) for
  masked similarity/distance computation. The shape of the matrix must
  be the same as the resulting matrix.

- min_simil:

  the minimum similarity value to be recorded.

- rank:

  an integer value specifying top-n most similarity values to be
  recorded.

- drop0:

  if `TRUE`, removes zero values to make the similarity/distance matrix
  sparse. It has no effect when `dense = TRUE`.

- diag:

  if `TRUE`, only compute diagonal elements of the similarity/distance
  matrix; useful when comparing corresponding rows or columns of `x` and
  `y`.

- use_nan:

  if `TRUE`, returns `NaN` if the standard deviation of a vector is zero
  when `method` is "correlation"; if all the values are zero in a vector
  when `method` is "cosine", "chisquared", "kullback", "jeffreys" or
  "jensen". Note that use of `NaN` makes the similarity/distance matrix
  denser and therefore larger in RAM. If `FALSE`, return zero in same
  use situations as above. If `NULL`, will also return zero but also
  generate a warning (default).

- sparse:

  if `TRUE`, returns
  [Matrix::sparseMatrix](https://rdrr.io/pkg/Matrix/man/sparseMatrix.html)
  object. When neither `min_simil` nor `rank` is used, dense matrices
  require less space in RAM.

- digits:

  determines rounding of small values towards zero. Use primarily to
  correct floating point errors. Rounding is performed in C++ in a
  similar way as [base::zapsmall](https://rdrr.io/r/base/zapsmall.html).

- p:

  weight for Minkowski distance.

- smooth:

  adds a fixed value to all the cells to avoid division by zero. Only
  used when `method` is "chisquared", "kullback", "jeffreys" or
  "jensen".

## Value

a [Matrix::Matrix](https://rdrr.io/pkg/Matrix/man/Matrix.html) object

## Details

### Available Methods

Similarity:

- `cosine`: cosine similarity

- `correlation`: Pearson's correlation

- `jaccard`: Jaccard coefficient

- `ejaccard`: the real value version of `jaccard`

- `fjaccard`: Fuzzy Jaccard coefficient

- `dice`: Dice coefficient

- `edice`: the real value version of `dice`

- `hamann`: Hamann similarity

- `faith`: Faith similarity

- `simple matching`: the percentage of common elements

Distance:

- `euclidean`: Euclidean distance

- `chisquared`: chi-squared distance

- `kullback`: Kullback–Leibler divergence

- `jeffreys`: Jeffreys divergence

- `jensen`: Jensen–Shannon divergence

- `manhattan`: Manhattan distance

- `maximum`: the largest difference between values

- `canberra`: Canberra distance

- `minkowski`: Minkowski distance

- `hamming`: Hamming distance

See the vignette for how the similarity and distance are computed:
[`vignette("measures", package = "proxyC")`](https://koheiw.github.io/proxyC/articles/measures.md)

### Parallel Computing

It performs parallel computing using Intel oneAPI Threads Building
Blocks. The number of threads for parallel computing should be specified
via `options(proxyC.threads)` before calling the functions. If the value
is -1, all the available threads will be used. Unless the option is
used, the number of threads will be limited by the environmental
variables (`OMP_THREAD_LIMIT` or `RCPP_PARALLEL_NUM_THREADS`) to comply
with CRAN policy and offer backward compatibility.

## See also

zapsmall

## Examples

``` r
mt <- Matrix::rsparsematrix(100, 100, 0.01)
simil(mt, method = "cosine")[1:5, 1:5]
#> Warning: x or y has vectors with all zero; consider setting use_nan = TRUE to set these values to NaN or use_nan = FALSE to suppress this warning
#> 5 x 5 sparse Matrix of class "dsTMatrix"
#>               
#> [1,] . . . . .
#> [2,] . 1 0 0 0
#> [3,] . 0 1 0 0
#> [4,] . 0 0 1 0
#> [5,] . 0 0 0 1
mt <- Matrix::rsparsematrix(100, 100, 0.01)
dist(mt, method = "euclidean")[1:5, 1:5]
#> 5 x 5 sparse Matrix of class "dsTMatrix"
#>                                                  
#> [1,] 0.000000 2.680019 2.830300 3.586712 3.074183
#> [2,] 2.680019 0.000000 0.910000 2.383695 1.506021
#> [3,] 2.830300 0.910000 0.000000 2.551490 1.759602
#> [4,] 3.586712 2.383695 2.551490 0.000000 2.819592
#> [5,] 3.074183 1.506021 1.759602 2.819592 0.000000
```
