# proxyC 0.2.0

## New features
- Added a `diag` argument to compute similarity/distance only for corresponding rows or columns (#13).
- Added a `smooth` parameter to chisquared and kullback leibler distances to solve negative values in sparse matrices (#15).
- Added the hamming distance (#18)

## Bug fixes
- Fixed the chi-squared distance to match `stats::chisq.test()` (#14).
- Fixed a bug in pairwise similarity/distance computation when `drop0 = TRUE` (#17).


# proxyC 0.1.5

## New features

- Add the `drop0` argument to address the floating point precision issue (#10).

## Bug fixes

- The digit argument is now passed to `dist()` (#11).


# proxyC 0.1.4

## New features

- Added `rowSds()`, `colSds()`, `rowZeros()` and `colZeros()` (#9).


# proxyC 0.1.3

## Bug fixes

- No longer assumes symmetry of resulting matrix when `x != y` (#4).

## New features

- Added the `digits` argument to correct rounding errors in C++ (#5).
