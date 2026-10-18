smt1 <- Matrix::rsparsematrix(100000, 100, 0.9)
smt2 <- Matrix::rsparsematrix(100000, 200, 0.9)
mt1 <- as.matrix(smt1)
mt2 <- as.matrix(smt2)
mask <- Matrix::rsparsematrix(100, 200, 1)

microbenchmark::microbenchmark(
    sparse = proxyC:::cpp_linear(smt1, smt2, 4, mask, 100, symm = TRUE, use_mask = TRUE, drop0 = TRUE),
    dense = proxyC:::cpp_linear(mt1, mt2, 4, mask, 100, symm = TRUE, use_mask = TRUE, drop0 = TRUE),
    times = 10
)

microbenchmark::microbenchmark(
    Matrix = Matrix::crossprod(smt1, smt2),
    proxyC = proxyC::crossprod(smt1, smt2),
    proxyC2 = proxyC::crossprod(smt1, smt2, sparse = FALSE),
    times = 10
)

microbenchmark::microbenchmark(
    Matrix = Matrix::crossprod(mt1, mt2),
    proxyC = proxyC::crossprod(mt1, mt2),
    proxyC2 = proxyC::crossprod(mt1, mt2, sparse = FALSE),
    times = 10
)

microbenchmark::microbenchmark(
    maxCol(smt1),
    max.col(smt1),
    times = 10
)

microbenchmark::microbenchmark(
    maxCol(mt1),
    max.col(mt1),
    times = 10
)

