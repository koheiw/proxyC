# Changelog

## proxyC 0.5.2

CRAN release: 2025-04-25

### Bug fixes

- Update tests to pass CRAN check with Intel MKL.

## proxyC 0.5.1

CRAN release: 2025-04-24

### Bug fixes

- Add workaround for a compilation error when clang 20 is used
  ([\#65](https://github.com/koheiw/proxyC/issues/65)).

## proxyC 0.5.0

CRAN release: 2025-04-22

### New features and improvements

- Enable [`simil()`](https://koheiw.github.io/proxyC/reference/simil.md)
  nad [`dist()`](https://koheiw.github.io/proxyC/reference/simil.md) to
  perform masked similarity/distance computation. If a pattern matrix is
  given via the `mask` argument, it computes scores for selected pairs
  of rows or columns. Pattern matrices can be created using newly added
  [`mask()`](https://koheiw.github.io/proxyC/reference/mask.md)
  function.
- Calculate “dice” and “edice” using Armadillo’s linear algebra
  functions in
  [`simil()`](https://koheiw.github.io/proxyC/reference/simil.md). As a
  result, computation of these scores became as fast as “cosine”,
  “correlation” and “euclidean”.
- Add
  [`crossprod()`](https://koheiw.github.io/proxyC/reference/crossprod.md)
  and
  [`tcrossprod()`](https://koheiw.github.io/proxyC/reference/crossprod.md)
  using the same infrastructure as
  [`simil()`](https://koheiw.github.io/proxyC/reference/simil.md) and
  [`dist()`](https://koheiw.github.io/proxyC/reference/simil.md).

## proxyC 0.4.2

### New features and improvements

- Reduce the overhead for dense similarity matrices by improving
  rounding numbers and conversion to Rcpp vectors.
- Return a dense matrix if `sparse = FALSE` to save space in RAM.

## proxyC 0.4.1

CRAN release: 2024-04-07

### Bug fixes

- Make detection of Intel oneAPI Threads Building Blocks (TBB) library
  more reliable.

## proxyC 0.4.0

CRAN release: 2024-04-05

### New features and improvements

- Use more recent Intel oneAPI Threads Building Blocks (TBB) library to
  improve the stability in parallel computing.
- Add `options(proxyC.threads)` to control the number of threads in
  parallel computing (but `RCPP_PARALLEL_NUM_THREADS` still has effect).

### New system requirements

- The RcppParallel package is no longer required as the TBB library in
  the operating system (Linux and MacOS) or Rtools (Windows) is used.
- Linux and MacOS must have the TBB library to enable parallel computing
  before installing this package from the source.

## proxyC 0.3.4

CRAN release: 2023-10-25

### New features and improvements

- Add “fjaccard” to
  [`simil()`](https://koheiw.github.io/proxyC/reference/simil.md) for
  Fuzzy Jaccard similarity
  ([\#42](https://github.com/koheiw/proxyC/issues/42)).

## proxyC 0.3.3

CRAN release: 2022-10-06

### New features and improvements

- Explicitly setting `use_nan = FALSE` will suppress warnings in
  [`simil()`](https://koheiw.github.io/proxyC/reference/simil.md) and
  [`dist()`](https://koheiw.github.io/proxyC/reference/simil.md).
- Add vignettes to explain how the similarity and distance measures are
  computed.

## proxyC 0.3.2

CRAN release: 2022-09-05

### Bug fixes

- Make further changes for Matrix v1.4-2.

## proxyC 0.3.1

CRAN release: 2022-08-22

### New features and improvements

- Add “jensen” to
  [`dist()`](https://koheiw.github.io/proxyC/reference/simil.md) for
  Jensen-Shannon divergence as a symmetric version of Kullback-Leibler
  divergence.
- Change how `x` and `y` are coerced to dgCMatrix for Matrix v1.4-2.

## proxyC 0.3.0

CRAN release: 2022-08-05

### New features and improvements

- Add “jeffreys” to
  [`dist()`](https://koheiw.github.io/proxyC/reference/simil.md) for
  Jeffreys divergence. It is a symmetric version of Kullback-Leibler
  divergence ([\#31](https://github.com/koheiw/proxyC/issues/31)).

## proxyC 0.2.4

CRAN release: 2021-12-10

### New features and improvements

- [`rowSds()`](https://koheiw.github.io/proxyC/reference/colSds.md),
  [`colSds()`](https://koheiw.github.io/proxyC/reference/colSds.md),
  [`rowZeros()`](https://koheiw.github.io/proxyC/reference/colZeros.md)
  and
  [`colZeros()`](https://koheiw.github.io/proxyC/reference/colZeros.md)
  return row or column names. They also work with both dense and sparse
  matrices ([\#28](https://github.com/koheiw/proxyC/issues/28)).

## proxyC 0.2.3

CRAN release: 2021-11-16

### Bug fixes

- Change “hamman” to “hamann” in
  [`simil()`](https://koheiw.github.io/proxyC/reference/simil.md) to
  correct misspelling
  ([\#26](https://github.com/koheiw/proxyC/issues/26)).

## proxyC 0.2.2

CRAN release: 2021-10-27

### New features

- [`simil()`](https://koheiw.github.io/proxyC/reference/simil.md) and
  [`dist()`](https://koheiw.github.io/proxyC/reference/simil.md) work
  with both dense and sparse matrices.
- `use_nan = TRUE` can be used not only for correlation but for all the
  distance and similarity measures.

## proxyC 0.2.1

CRAN release: 2021-09-02

### New features

- Computing the correlation similarity on vectors with a standard
  deviation will generate a zero correlation and a warning. The warning
  can be turned off by setting `use_nan = TRUE`, in which case the
  computed correlation similarity will be `NaN` instead
  ([\#21](https://github.com/koheiw/proxyC/issues/21)).

### Bug fixes

- Fixed infinite values being returned by the correlation similarity
  ([\#20](https://github.com/koheiw/proxyC/issues/20)).

## proxyC 0.2.0

CRAN release: 2021-05-11

### New features

- Added a `diag` argument to compute similarity/distance only for
  corresponding rows or columns
  ([\#13](https://github.com/koheiw/proxyC/issues/13)).
- Added a `smooth` parameter to chisquared and kullback leibler
  distances to solve negative values in sparse matrices
  ([\#15](https://github.com/koheiw/proxyC/issues/15)).
- Added the hamming distance
  ([\#18](https://github.com/koheiw/proxyC/issues/18))

### Bug fixes

- Fixed the chi-squared distance to match
  [`stats::chisq.test()`](https://rdrr.io/r/stats/chisq.test.html)
  ([\#14](https://github.com/koheiw/proxyC/issues/14)).
- Fixed a bug in pairwise similarity/distance computation when
  `drop0 = TRUE` ([\#17](https://github.com/koheiw/proxyC/issues/17)).

## proxyC 0.1.5

CRAN release: 2019-07-21

### New features

- Add the `drop0` argument to address the floating point precision issue
  ([\#10](https://github.com/koheiw/proxyC/issues/10)).

### Bug fixes

- The digit argument is now passed to
  [`dist()`](https://koheiw.github.io/proxyC/reference/simil.md)
  ([\#11](https://github.com/koheiw/proxyC/issues/11)).

## proxyC 0.1.4

CRAN release: 2019-06-04

### New features

- Add [`rowSds()`](https://koheiw.github.io/proxyC/reference/colSds.md),
  [`colSds()`](https://koheiw.github.io/proxyC/reference/colSds.md),
  [`rowZeros()`](https://koheiw.github.io/proxyC/reference/colZeros.md)
  and
  [`colZeros()`](https://koheiw.github.io/proxyC/reference/colZeros.md)
  ([\#9](https://github.com/koheiw/proxyC/issues/9)).

## proxyC 0.1.3

CRAN release: 2019-04-21

### Bug fixes

- No longer assumes symmetry of resulting matrix when `x != y`
  ([\#4](https://github.com/koheiw/proxyC/issues/4)).

### New features

- Add the `digits` argument to correct rounding errors in C++
  ([\#5](https://github.com/koheiw/proxyC/issues/5)).
