require(Matrix)
source("../../tests/function.R")

mat_test <- rsparsematrix(100, 50, 0.5)
mat_test[1, ] <- 0.5 # add one row with sd(x) == 0
mat_test[, 1] <- 0.5 # add one col with sd(x) == 0
mat_test[2, ] <- 0.0 # add one row with sum(x) == 0
mat_test[, 2] <- 0.0 # add one col with sum(x) == 0

test_that("test cosine similarity", {
    skip_if_not_installed("proxy")
    test_simil(mat_test, "cosine", margin = 1, use_nan = TRUE)
    test_simil(mat_test, "cosine", margin = 2, use_nan = TRUE)
})

test_that("test correlation similarity", {
    skip_if_not_installed("proxy")
    test_simil(mat_test, "correlation", margin = 1, use_nan = TRUE)
    test_simil(mat_test, "correlation", margin = 2, use_nan = TRUE)
})

test_that("test jaccard similarity", {
    skip_if_not_installed("proxy")
    # proxy::simil(method = "jaccard") wrongly returns 1 for all zero vector
    test_simil(mat_test[-2, -2], "jaccard", margin = 1)
    test_simil(mat_test[-2, -2], "jaccard", margin = 2)
})

test_that("test ejaccard similarity", {
    skip_if_not_installed("proxy")
    # proxy::simil(method = "ejaccard") wrongly returns 1 for all0  vector
    test_simil(mat_test[-2, -2], "ejaccard", margin = 1)
    test_simil(mat_test[-2, -2], "ejaccard", margin = 2)
})

test_that("test dice similarity", {
    skip_if_not_installed("proxy")
    # proxy::simil(method = "dice") wrongly returns 1 for all0 or sd0 vector
    test_simil(mat_test[c(-1, -2), c(-1, -2)], "dice", margin = 1)
    test_simil(mat_test[c(-1, -2), c(-1, -2)], "dice", margin = 2)
})

test_that("test edice similarity", {
    skip_if_not_installed("proxy")
    # proxy::simil(method = "edice") wrongly returns 1 for all0 or sd0 vector
    test_simil(mat_test[c(-1, -2), c(-1, -2)], "edice", margin = 1)
    test_simil(mat_test[c(-1, -2), c(-1, -2)], "edice", margin = 2)
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
    expect_warning(proxyC::simil(mat1, mat4, method = "cosine", use_nan = FALSE),
                   "x or y has vectors with all zero; consider setting use_nan = TRUE")
    suppressWarnings({
        expect_equal(proxyC::simil(mat1, mat2, method = "correlation", use_nan = FALSE)[1,1], 0)
        expect_equal(proxyC::simil(mat1, mat3, method = "correlation", use_nan = FALSE)[1,1], 0)
        expect_equal(proxyC::simil(mat1, mat2, method = "correlation", use_nan = FALSE, diag = TRUE)[1,1], 0)
        expect_equal(proxyC::simil(mat1, mat3, method = "correlation", use_nan = FALSE, diag = TRUE)[1,1], 0)
    })

    expect_true(is.na(proxyC::simil(mat1, mat2, method = "correlation", use_nan = TRUE)[1,1]))
    expect_true(is.na(proxyC::simil(mat1, mat3, method = "correlation", use_nan = TRUE)[1,1]))
    expect_true(is.na(proxyC::simil(mat1, mat2, method = "correlation", use_nan = TRUE, diag = TRUE)[1,1]))
    expect_true(is.na(proxyC::simil(mat1, mat3, method = "correlation", use_nan = TRUE, diag = TRUE)[1,1]))

    expect_false(is.nan(proxyC::simil(mat1, mat2, method = "cosine", use_nan = TRUE)[1,1]))
    expect_false(is.nan(proxyC::simil(mat1, mat3, method = "cosine", use_nan = TRUE)[1,1]))
    expect_true(is.nan(proxyC::simil(mat1, mat4, method = "cosine", use_nan = TRUE)[1,1]))
})


test_that("simil returns zero or NaN correctly", {

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

})
