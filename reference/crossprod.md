# Cross-product of large sparse matrices

Compute the (transposed) cross-product of large sparse matrices using
the same infrastructure as
[`simil()`](https://koheiw.github.io/proxyC/reference/simil.md) and
[`dist()`](https://koheiw.github.io/proxyC/reference/simil.md).

## Usage

``` r
crossprod(x, y = NULL, min_prod = NULL, sparse = TRUE, digits = 14)

tcrossprod(x, y = NULL, min_prod = NULL, sparse = TRUE, digits = 14)
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

- min_prod:

  the minimum product to be recorded.

- sparse:

  if `TRUE`, returns
  [Matrix::sparseMatrix](https://rdrr.io/pkg/Matrix/man/sparseMatrix.html)
  object. When neither `min_simil` nor `rank` is used, dense matrices
  require less space in RAM.

- digits:

  determines rounding of small values towards zero. Use primarily to
  correct floating point errors. Rounding is performed in C++ in a
  similar way as [base::zapsmall](https://rdrr.io/r/base/zapsmall.html).

## Value

a [Matrix::Matrix](https://rdrr.io/pkg/Matrix/man/Matrix.html) object
