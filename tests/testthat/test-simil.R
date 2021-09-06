require(Matrix)
mat_test <- rsparsematrix(100, 50, 0.5)
mat_test[1, ] <- 0.5 # add one row with sd(x) == 0
mat_test[, 1] <- 0.5 # add one col with sd(x) == 0
mat_test[2, ] <- 0.0 # add one row with sum(x) == 0
mat_test[, 2] <- 0.0 # add one col with sum(x) == 0

# returns TRUE where values are affected by vectors with all zero
is_all0 <- function(x, y = x, margin = 1) {
    if (margin == 1) {
        r1 <- outer(rowZeros(x) == ncol(x), rep(TRUE, nrow(y)))
        r2 <- outer(rep(TRUE, nrow(x)), rowZeros(y) == ncol(y))
    } else {
        r1 <- outer(colZeros(x) == nrow(x), rep(TRUE, ncol(y)))
        r2 <- outer(rep(TRUE, ncol(x)), colZeros(y) == nrow(y))
    }
    return(r1 | r2)
}

# returns TRUE where values are affected by vectors with zero variance
is_sd0 <- function(x, y = x, margin = 1) {
    if (margin == 1) {
        r1 <- outer(rowSds(x) == 0, rep(TRUE, nrow(y)))
        r2 <- outer(rep(TRUE, nrow(x)), rowSds(y) == 0)
    } else {
        r1 <- outer(colSds(x) == 0, rep(TRUE, ncol(y)))
        r2 <- outer(rep(TRUE, ncol(x)), colSds(y) == 0)
    }
    return(r1 | r2)
}

test_simil <- function(x, method, margin, ignore_upper = FALSE, ignore_diag = TRUE, use_nan = FALSE, ...) {
    # test with only x
    suppressWarnings({
    s1 <- as.matrix(simil(x, method = method, margin = margin, use_nan = use_nan, ...))
    })
    s2 <- proxy::as.matrix(proxy::simil(as.matrix(x),
                           method = method, by_rows = margin == 1, diag = TRUE, ...))
    if (ignore_diag)
        diag(s1) <- diag(s2) <- 0
    if (ignore_upper)
        s1[upper.tri(s1, TRUE)] <- s2[upper.tri(s2, TRUE)] <- 0
    if (use_nan) {
        s2[is_all0(x, y, margin = margin)] <- NaN
        if (method == "correlation")
            s2[is_sd0(x, y, margin = margin)] <- NaN
    }
    expect_equal(as.numeric(s1), as.numeric(s2), tolerance = 0.001)

    # test with x and y, different size
    if (margin == 1) {
        y <- x[sample(nrow(x), 10),]
    } else {
        y <- x[,sample(ncol(x), 10)]
    }
    suppressWarnings({
    s3 <- as.matrix(simil(x, y, method = method, margin = margin, use_nan = use_nan, ...))
    })
    s4 <- as.matrix(proxy::simil(as.matrix(x), as.matrix(y),
                                 method = method, by_rows = margin == 1, diag = TRUE, ...))
    if (ignore_diag)
        diag(s3) <- diag(s4) <- 0
    if (ignore_upper)
        s3[upper.tri(s3, TRUE)] <- s4[upper.tri(s4, TRUE)] <- 0
    if (use_nan) {
        s3[is_all0(x, y, margin = margin)] <- NaN
        if (method == "correlation")
            s3[is_sd0(x, y, margin = margin)] <- NaN
    }
    expect_equal(as.numeric(s3), as.numeric(s4), tolerance = 0.001)


    # test with x and y, same size
    if (margin == 1) {
        y <- x[sample(nrow(x)),]
    } else {
        y <- x[,sample(ncol(x))]
    }

    suppressWarnings({
    s5 <- as.matrix(simil(x, y, method = method, margin = margin, use_nan = use_nan, ...))
    })
    s6 <- as.matrix(proxy::simil(as.matrix(x), as.matrix(y),
                                 method = method, by_rows = margin == 1, diag = TRUE, ...))
    if (ignore_diag)
        diag(s5) <- diag(s6) <- 0
    if (ignore_upper)
        s5[upper.tri(s5, TRUE)] <- s6[upper.tri(s6, TRUE)] <- 0
    if (use_nan) {
        s6[is_all0(x, y, margin = margin)] <- NaN
        if (method == "correlation")
            s6[is_sd0(x, y, margin = margin)] <- NaN
    }
    print(s5)
    print(s6)
    expect_equal(as.numeric(s5), as.numeric(s6), tolerance = 0.001)
}

test_that("test cosine similarity", {
    skip_if_not_installed("proxy")
    test_simil(mat_test, "cosine", margin = 1)
    test_simil(mat_test, "cosine", margin = 2)
})

test_that("test correlation similarity", {
    skip_if_not_installed("proxy")
    test_simil(mat_test, "correlation", margin = 1, use_nan = TRUE)
    test_simil(mat_test, "correlation", margin = 2, use_nan = TRUE)
})

test_that("test jaccard similarity", {
    skip_if_not_installed("proxy")
    test_simil(mat_test, "jaccard", margin = 1)
    test_simil(mat_test, "jaccard", margin = 2)
})

test_that("test ejaccard similarity", {
    skip_if_not_installed("proxy")
    test_simil(mat_test, "ejaccard", margin = 1)
    test_simil(mat_test, "ejaccard", margin = 2)
})

test_that("test dice similarity", {
    skip_if_not_installed("proxy")
    test_simil(mat_test, "dice", margin = 1)
    test_simil(mat_test, "dice", margin = 2)
})

test_that("test edice similarity", {
    skip_if_not_installed("proxy")
    test_simil(mat_test, "edice", margin = 1)
    test_simil(mat_test, "edice", margin = 2)
})

test_that("test simple matching similarity", {
    skip_if_not_installed("proxy")
    test_simil(mat_test, "simple matching", margin = 1)
    test_simil(mat_test, "simple matching", margin = 2)
})

test_that("test hamman similarity", {
    skip_if_not_installed("proxy")
    test_simil(mat_test, "hamman", margin = 1)
    test_simil(mat_test, "hamman", margin = 2)
})

test_that("test faith similarity", {
    skip_if_not_installed("proxy")
    test_simil(mat_test, "faith", margin = 1, ignore_upper = TRUE)
    test_simil(mat_test, "faith", margin = 2, ignore_upper = TRUE)
})


test_that("use_na is working", {

    mat1 <- Matrix::Matrix(1:4, nrow = 1, sparse = TRUE)
    mat2 <- Matrix::Matrix(rep(0.001, 4), nrow = 1, sparse = TRUE)
    mat3 <- Matrix::Matrix(rep(1, 4), nrow = 1, sparse = TRUE)
    mat4 <- Matrix::Matrix(rep(0, 4), nrow = 1, sparse = TRUE)

    expect_warning(proxyC::simil(mat1, mat2, method = "correlation", use_nan = FALSE),
                   "x or y has vectors with zero standard deviation; consider setting use_nan = TRUE")
    suppressWarnings({
        expect_equal(proxyC::simil(mat1, mat2, method = "correlation", use_nan = FALSE)[1,1], 0)
        expect_equal(proxyC::simil(mat1, mat3, method = "correlation", use_nan = FALSE)[1,1], 0)
        expect_equal(proxyC::simil(mat1, mat2, method = "correlation", use_nan = FALSE, diag = TRUE)[1,1], 0)
        expect_equal(proxyC::simil(mat1, mat3, method = "correlation", use_nan = FALSE, diag = TRUE)[1,1], 0)
    })

    expect_equal(proxyC::simil(mat1, mat2, method = "correlation", use_nan = TRUE)[1,1], NaN)
    expect_equal(proxyC::simil(mat1, mat3, method = "correlation", use_nan = TRUE)[1,1], NaN)
    expect_equal(proxyC::simil(mat1, mat2, method = "correlation", use_nan = TRUE, diag = TRUE)[1,1], NaN)
    expect_equal(proxyC::simil(mat1, mat3, method = "correlation", use_nan = TRUE, diag = TRUE)[1,1], NaN)

    expect_false(is.nan(proxyC::simil(mat1, mat2, method = "cosine", use_nan = TRUE)[1,1]))
    expect_false(is.nan(proxyC::simil(mat1, mat3, method = "cosine", use_nan = TRUE)[1,1]))
    expect_false(is.nan(proxyC::simil(mat1, mat4, method = "cosine", use_nan = TRUE)[1,1]))

})


# TODO: test there is no zero or inf

mat <- Matrix::Matrix(matrix(c(0, 0, 0,
                               1, 1, 1,
                               1, 5, 2,
                               2, 3, 4), byrow = TRUE, nrow = 4), sparse = TRUE)

expect_equivalent(
    as.matrix(proxyC::simil(mat, metho = "cosine", margin = 1) == 0),
    is_all0(mat, margin = 1)
)

expect_equivalent(
    as.matrix(proxyC::simil(mat, metho = "cosine", margin = 2) == 0),
    is_all0(mat, margin = 2)
)

expect_equivalent(
    suppressWarnings(as.matrix(proxyC::simil(mat, method = "correlation", margin = 1) == 0)),
    is_all0(mat, margin = 1) | is_sd0(mat, margin = 1)
)

expect_equivalent(
    suppressWarnings(as.matrix(proxyC::simil(mat, method = "correlation", margin = 2) == 0)),
    is_all0(mat, margin = 2) | is_sd0(mat, margin = 2)
)

expect_equivalent(
    is.nan(as.matrix(proxyC::simil(mat, metho = "cosine", margin = 1, use_nan = TRUE))),
    is_all0(mat, margin = 1)
)

expect_equivalent(
    is.nan(as.matrix(proxyC::simil(mat, metho = "cosine", margin = 2, use_nan = TRUE))),
    is_all0(mat, margin = 2)
)

expect_equivalent(
    is.nan(as.matrix(proxyC::simil(mat, method = "correlation", margin = 1, use_nan = TRUE))),
    is_all0(mat, margin = 1) | is_sd0(mat, margin = 1)
)

expect_equivalent(
    is.nan(as.matrix(proxyC::simil(mat, method = "correlation", margin = 2, use_nan = TRUE))),
    is_all0(mat, margin = 2) | is_sd0(mat, margin = 2)
)

is_all0(mat, margin = 1)
is_all0(mat, margin = 2)
is_sd0(mat, margin = 1)
is_sd0(mat, margin = 2)

is_all0(mat, mat[1:3,], margin = 1)
is_sd0(mat, mat[1:3,], margin = 1)
mat
mat[1:3,]

is_all0(mat, mat[,1:2], margin = 2)
is_sd0(mat, mat[,1:2], margin = 2)

proxyC::simil(mat, method = "cosine")
proxyC::simil(mat, method = "correlation")
proxyC::simil(mat, method = "jaccard")

proxyC::simil(mat, method = "jaccard", use_nan = TRUE)



proxyC::simil(mat, method = "correlation")
cor(t(as.matrix(mat)))


proxyC::simil(mat, method = "correlation")
proxyC::simil(mat, method = "correlation", diag = TRUE)
proxyC::simil(mat, method = "cosine", diag = TRUE)
proxyC::dist(mat, method = "euclidean", diag = TRUE)

proxyC::simil(mat, method = "jaccard")
proxyC::simil(mat, method = "jaccard", use_nan = TRUE)


proxyC::simil(mat, method = "dice")
proxyC::simil(mat, method = "hamman")
proxyC::simil(mat, method = "simple matching")
proxyC::simil(mat, method = "faith")

cor(t(as.matrix(mat)))
