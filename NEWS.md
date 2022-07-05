# proxyC 0.3.0

## New features and improvements

- Add "jeffreys" to `dist()` for Jeffreys divergence. It is a symmetric version of Kullback-Leibluer divergence (#31).

# proxyC 0.2.4

## New features and improvements

- `rowSds()`, `colSds()`, `rowZeros()` and `colZeros()` return row or column names. They also work with both dense and sparse matrices (#28).

# proxyC 0.2.3

## Bug fixes

- Change "hamman" to "hamann" in `simil()` to correct misspelling (#26).

# proxyC 0.2.2

## New features

- `simil()` and `dist()` work with both dense and sparse matrices.
- `use_nan = TRUE` can be used not only for correlation but for all the distance 
  and similarity measures.

# proxyC 0.2.1

## New features

- Computing the correlation similarity on vectors with a standard deviation will 
  generate a zero correlation and a warning. The warning can be turned off by 
  setting `use_nan = TRUE`, in which case the computed correlation similarity 
  will be `NaN` instead (#21).

## Bug fixes

- Fixed infinite values being returned by the correlation similarity (#20).

# proxyC 0.2.0

## New features

- Added a `diag` argument to compute similarity/distance only for corresponding 
  rows or columns (#13).
- Added a `smooth` parameter to chisquared and kullback leibler distances to 
  solve negative values in sparse matrices (#15).
- Added the hamming distance (#18)

## Bug fixes

- Fixed the chi-squared distance to match `stats::chisq.test()` (#14).
- Fixed a bug in pairwise similarity/distance computation when `drop0 = TRUE` 
  (#17).

# proxyC 0.1.5

## New features

- Add the `drop0` argument to address the floating point precision issue (#10).

## Bug fixes

- The digit argument is now passed to `dist()` (#11).

# proxyC 0.1.4

## New features

- Add `rowSds()`, `colSds()`, `rowZeros()` and `colZeros()` (#9).

# proxyC 0.1.3

## Bug fixes

- No longer assumes symmetry of resulting matrix when `x != y` (#4).

## New features

- Add the `digits` argument to correct rounding errors in C++ (#5).
