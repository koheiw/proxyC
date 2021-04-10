# v0.2.0

- Add the `diag` argument to compute similarity/distance only for corresponding rows or columns.
- Correct the chi-squared distance to match `stats::chisq.test`.

# v0.1.5

## New feature

- Add the dorp0 argument to address the floating point precision issue

## Bug fix

- The digit argument is now passed to `dist()`.

# v0.1.4

## New feature

- Add `rowSds()`, `colSds()`, `rowZeros()` and `colZeros()`.

# v0.1.3

## Bug fix

- No longer assumes symmetry of resulting matrix when `x != y`.

## New feature

- Add the digits argument to correct rounding errors in C++ (#5).
