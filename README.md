
# proxyC: R package for large-scale similarity/distance computation

**proxyC** computes proximity between rows or columns of large matrices
efficiently in C++. It is optimized for large sparse matrices using the
Armadillo and Intel TBB libraries. Among several built-in
similarity/distance measures, computation of correlation, cosine
similarity and Euclidean distance is particularly fast.

This code was originally written for
[**quanteda**](https://github.com/quanteda/quanteda) to compute
similarity/distance between documents or features in large corpora, but
separated as a stand-alone package to make it available for broader data
scientific purposes.

``` r
install.packages("proxyC")
```

``` r
require(Matrix)
## Loading required package: Matrix
require(microbenchmark)
## Loading required package: microbenchmark
require(RcppParallel)
## Loading required package: RcppParallel
require(ggplot2)
## Loading required package: ggplot2

# Set number of threads
setThreadOptions(8)

# Make a matrix with 99% zeros
sm1k <- rsparsematrix(1000, 1000, 0.01) # 1,000 columns
sm10k <- rsparsematrix(1000, 10000, 0.01) # 10,000 columns

# Convert to dense format
dm1k <- as.matrix(sm1k) 
dm10k <- as.matrix(sm10k)
```

## Cosine similarity between columns

With sparse matrices, **proxyC** is roughly 10 to 100 times faster than
**proxy**.

``` r
bm1 <- microbenchmark(
    "proxyC 1k" = proxyC::simil(sm1k, margin = 2, method = "cosine"),
    "proxy 1k" = proxy::simil(dm1k, method = "cosine"),
    "proxyC 10k" = proxyC::simil(sm10k, margin = 2, method = "cosine"),
    "proxy 10k" = proxy::simil(dm10k, method = "cosine"),
    times = 10
)
autoplot(bm1)
## Coordinate system already present. Adding new coordinate system, which will replace the existing one.
```

![](man/images/unnamed-chunk-4-1.png)<!-- -->

## Top-10 cosine similarity

If `rank` is used, **proxyC** becomes even faster as many similarity
scores are discarded (rounded to zero).

``` r
bm2 <- microbenchmark(
    "proxyC rank" = proxyC::simil(sm1k, margin = 2, method = "cosine", rank = 10),
    "proxyC all" = proxyC::simil(sm1k, margin = 2, method = "cosine"),
    times = 10
)
autoplot(bm2)
## Coordinate system already present. Adding new coordinate system, which will replace the existing one.
```

![](man/images/unnamed-chunk-5-1.png)<!-- -->

## Correlation greater than 0.9

`min_simil` also makes **proxyC** faster.

``` r
bm3 <- microbenchmark(
    "proxyC min_simil" = proxyC::simil(sm1k, margin = 2, method = "correlation", min_simil = 0.9),
    "proxyC all" = proxyC::simil(sm1k, margin = 2, method = "correlation"),
    times = 10
)
autoplot(bm3)
## Coordinate system already present. Adding new coordinate system, which will replace the existing one.
```

![](man/images/unnamed-chunk-6-1.png)<!-- -->
