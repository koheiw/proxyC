source("function.R")

test_that("test test objects and functions", {
    expect_equal(
        rowSds(mat_test)[1:3],
        c(0, 0, 0)
    )
    expect_equal(
        rowZeros(mat_test)[1:3],
        c(50, 0, 0)
    )

    mat_sd0 <- mat_all0 <- matrix(FALSE, nrow(mat_test), nrow(mat_test))
    mat_all0[1,] <- mat_all0[,1]<- TRUE
    expect_equivalent(
        is_all0(mat_test), mat_all0
    )
    mat_sd0[1:3,]  <- mat_sd0[,1:3] <- TRUE
    expect_equivalent(
        is_sd0(mat_test), mat_sd0
    )
})

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
    test_simil(mat_test[-1,], "jaccard", margin = 1)
    test_simil(mat_test[-1,], "jaccard", margin = 2)
})

test_that("test ejaccard similarity", {
    skip_if_not_installed("proxy")
    # proxy::simil(method = "ejaccard") wrongly returns 1 for all0  vector
    test_simil(mat_test[-1,], "ejaccard", margin = 1)
    test_simil(mat_test[-1,], "ejaccard", margin = 2)
})

test_that("test dice similarity", {
    skip_if_not_installed("proxy")
    # proxy::simil(method = "dice") wrongly returns 1 for all0 or sd0 vector
    test_simil(mat_test[1:3 * -1,], "dice", margin = 1)
    test_simil(mat_test[1:3 * -1,], "dice", margin = 2)
})

test_that("test edice similarity", {
    skip_if_not_installed("proxy")
    # proxy::simil(method = "edice") wrongly returns 1 for all0 or sd0 vector
    test_simil(mat_test[1:3 * -1,], "edice", margin = 1)
    test_simil(mat_test[1:3 * -1,], "edice", margin = 2)
})

test_that("test simple matching similarity", {
    skip_if_not_installed("proxy")
    test_simil(mat_test, "simple matching", margin = 1)
    test_simil(mat_test, "simple matching", margin = 2)
})

test_that("test hamanm similarity", {
    skip_if_not_installed("proxy")
    test_simil(mat_test, "hamman", margin = 1) # use misspelling for proxy
    test_simil(mat_test, "hamman", margin = 2) # use misspelling for proxy
})

test_that("test faith similarity", {
    skip_if_not_installed("proxy")
    test_simil(mat_test, "faith", margin = 1, ignore_upper = TRUE)
    test_simil(mat_test, "faith", margin = 2, ignore_upper = TRUE)
})


test_that("use_nan is working", {

    mat1 <- Matrix::Matrix(1:4, nrow = 1, sparse = TRUE)
    mat2 <- Matrix::Matrix(rep(0.0001, 4), nrow = 1, sparse = TRUE)
    mat3 <- Matrix::Matrix(rep(1, 4), nrow = 1, sparse = TRUE)
    mat4 <- Matrix::Matrix(rep(0, 4), nrow = 1, sparse = TRUE)

    expect_warning(proxyC::simil(mat1, mat2, method = "correlation", use_nan = FALSE),
                   "x or y has vectors with zero standard deviation; consider setting use_nan = TRUE")
    expect_warning(proxyC::simil(mat1, mat4, method = "cosine", use_nan = FALSE),
                   "x or y has vectors with all zero; consider setting use_nan = TRUE")
    suppressWarnings({
        expect_equal(proxyC::simil(mat1, mat2, method = "correlation", use_nan = FALSE)[1,1], 0)
        expect_equal(proxyC::simil(mat1, mat3, method = "correlation", use_nan = FALSE)[1,1], 0)
        expect_equal(proxyC::simil(mat1, mat4, method = "correlation", use_nan = FALSE)[1,1], 0)
        expect_equal(proxyC::simil(mat2, mat1, method = "correlation", use_nan = FALSE)[1,1], 0)
        expect_equal(proxyC::simil(mat3, mat1, method = "correlation", use_nan = FALSE)[1,1], 0)
        expect_equal(proxyC::simil(mat4, mat1, method = "correlation", use_nan = FALSE)[1,1], 0)

        expect_equal(proxyC::simil(mat1, mat2, method = "correlation", use_nan = FALSE, diag = TRUE)[1,1], 0)
        expect_equal(proxyC::simil(mat1, mat3, method = "correlation", use_nan = FALSE, diag = TRUE)[1,1], 0)
        expect_equal(proxyC::simil(mat1, mat4, method = "correlation", use_nan = FALSE, diag = TRUE)[1,1], 0)
        expect_equal(proxyC::simil(mat2, mat1, method = "correlation", use_nan = FALSE, diag = TRUE)[1,1], 0)
        expect_equal(proxyC::simil(mat3, mat1, method = "correlation", use_nan = FALSE, diag = TRUE)[1,1], 0)
        expect_equal(proxyC::simil(mat4, mat1, method = "correlation", use_nan = FALSE, diag = TRUE)[1,1], 0)
    })

    expect_true(is.nan(proxyC::simil(mat1, mat2, method = "correlation", use_nan = TRUE)[1,1]))
    expect_true(is.nan(proxyC::simil(mat1, mat3, method = "correlation", use_nan = TRUE)[1,1]))
    expect_true(is.nan(proxyC::simil(mat1, mat4, method = "correlation", use_nan = TRUE)[1,1]))
    expect_true(is.nan(proxyC::simil(mat2, mat1, method = "correlation", use_nan = TRUE)[1,1]))
    expect_true(is.nan(proxyC::simil(mat3, mat1, method = "correlation", use_nan = TRUE)[1,1]))
    expect_true(is.nan(proxyC::simil(mat4, mat1, method = "correlation", use_nan = TRUE)[1,1]))

    expect_true(is.nan(proxyC::simil(mat1, mat2, method = "correlation", use_nan = TRUE, diag = TRUE)[1,1]))
    expect_true(is.nan(proxyC::simil(mat1, mat3, method = "correlation", use_nan = TRUE, diag = TRUE)[1,1]))
    expect_true(is.nan(proxyC::simil(mat1, mat4, method = "correlation", use_nan = TRUE, diag = TRUE)[1,1]))
    expect_true(is.nan(proxyC::simil(mat2, mat1, method = "correlation", use_nan = TRUE, diag = TRUE)[1,1]))
    expect_true(is.nan(proxyC::simil(mat3, mat1, method = "correlation", use_nan = TRUE, diag = TRUE)[1,1]))
    expect_true(is.nan(proxyC::simil(mat4, mat1, method = "correlation", use_nan = TRUE, diag = TRUE)[1,1]))

    suppressWarnings({
        expect_equal(proxyC::simil(mat1, mat4, method = "cosine", use_nan = FALSE)[1,1], 0)
        expect_equal(proxyC::simil(mat4, mat1, method = "cosine", use_nan = FALSE)[1,1], 0)

        expect_equal(proxyC::simil(mat1, mat4, method = "cosine", use_nan = FALSE, diag = TRUE)[1,1], 0)
        expect_equal(proxyC::simil(mat4, mat1, method = "cosine", use_nan = FALSE, diag = TRUE)[1,1], 0)
    })
    expect_false(is.nan(proxyC::simil(mat1, mat2, method = "cosine", use_nan = TRUE)[1,1]))
    expect_false(is.nan(proxyC::simil(mat1, mat3, method = "cosine", use_nan = TRUE)[1,1]))
    expect_true(is.nan(proxyC::simil(mat1, mat4, method = "cosine", use_nan = TRUE)[1,1]))
    expect_false(is.nan(proxyC::simil(mat2, mat1, method = "cosine", use_nan = TRUE)[1,1]))
    expect_false(is.nan(proxyC::simil(mat3, mat1, method = "cosine", use_nan = TRUE)[1,1]))
    expect_true(is.nan(proxyC::simil(mat4, mat1, method = "cosine", use_nan = TRUE)[1,1]))

    expect_false(is.nan(proxyC::simil(mat1, mat2, method = "cosine", use_nan = TRUE, diag = TRUE)[1,1]))
    expect_false(is.nan(proxyC::simil(mat1, mat3, method = "cosine", use_nan = TRUE, diag = TRUE)[1,1]))
    expect_true(is.nan(proxyC::simil(mat1, mat4, method = "cosine", use_nan = TRUE, diag = TRUE)[1,1]))
    expect_false(is.nan(proxyC::simil(mat2, mat1, method = "cosine", use_nan = TRUE, diag = TRUE)[1,1]))
    expect_false(is.nan(proxyC::simil(mat3, mat1, method = "cosine", use_nan = TRUE, diag = TRUE)[1,1]))
    expect_true(is.nan(proxyC::simil(mat4, mat1, method = "cosine", use_nan = TRUE, diag = TRUE)[1,1]))

})


test_that("simil returns zero or NaN correctly", {

    mat <- Matrix::Matrix(matrix(c(0, 0, 0,
                                   1, 1, 1,
                                   1, 5, 2,
                                   2, 3, 4), byrow = TRUE, nrow = 4), sparse = TRUE)

    # cosine
    expect_equivalent(
        suppressWarnings(as.matrix(proxyC::simil(mat, method = "cosine", margin = 1) == 0)),
        is_all0(mat, margin = 1)
    )
    expect_equivalent(
        as.matrix(proxyC::simil(mat, method = "cosine", margin = 2) == 0),
        is_all0(mat, margin = 2)
    )
    expect_equivalent(
        is.nan(as.matrix(proxyC::simil(mat, method = "cosine", margin = 1, use_nan = TRUE))),
        is_all0(mat, margin = 1)
    )
    expect_equivalent(
        is.nan(as.matrix(proxyC::simil(mat, method = "cosine", margin = 2, use_nan = TRUE))),
        is_all0(mat, margin = 2)
    )

    # correlation
    expect_equivalent(
        suppressWarnings(as.matrix(proxyC::simil(mat, method = "correlation", margin = 1) == 0)),
        is_all0(mat, margin = 1) | is_sd0(mat, margin = 1)
    )
    expect_equivalent(
        as.matrix(proxyC::simil(mat, method = "correlation", margin = 2) == 0),
        is_all0(mat, margin = 2) | is_sd0(mat, margin = 2)
    )
    expect_equivalent(
        is.nan(as.matrix(proxyC::simil(mat, method = "correlation", margin = 1, use_nan = TRUE))),
        is_all0(mat, margin = 1) | is_sd0(mat, margin = 1)
    )
    expect_equivalent(
        is.nan(as.matrix(proxyC::simil(mat, method = "correlation", margin = 2, use_nan = TRUE))),
        is_all0(mat, margin = 2) | is_sd0(mat, margin = 2)
    )

    # jaccard
    expect_false(
        any(is.nan(as.numeric(proxyC::simil(mat, method = "jaccard", margin = 1, use_nan = TRUE))))
    )
    expect_false(
        any(is.nan(as.numeric(proxyC::simil(mat, method = "jaccard", margin = 2, use_nan = TRUE))))
    )

    # ejaccard
    expect_false(
        any(is.nan(as.numeric(proxyC::simil(mat, method = "ejaccard", margin = 1, use_nan = TRUE))))
    )
    expect_false(
        any(is.nan(as.numeric(proxyC::simil(mat, method = "ejaccard", margin = 2, use_nan = TRUE))))
    )

    # dice
    expect_false(
        any(is.nan(as.numeric(proxyC::simil(mat, method = "dice", margin = 1, use_nan = TRUE))))
    )
    expect_false(
        any(is.nan(as.numeric(proxyC::simil(mat, method = "dice", margin = 2, use_nan = TRUE))))
    )

    # edice
    expect_false(
        any(is.nan(as.numeric(proxyC::simil(mat, method = "edice", margin = 1, use_nan = TRUE))))
    )
    expect_false(
        any(is.nan(as.numeric(proxyC::simil(mat, method = "edice", margin = 2, use_nan = TRUE))))
    )

    # hamman
    expect_false(
        any(is.nan(as.numeric(proxyC::simil(mat, method = "hamann", margin = 1, use_nan = TRUE))))
    )
    expect_false(
        any(is.nan(as.numeric(proxyC::simil(mat, method = "hamann", margin = 2, use_nan = TRUE))))
    )

    # simple matching
    expect_false(
        any(is.nan(as.numeric(proxyC::simil(mat, method = "simple matching", margin = 1, use_nan = TRUE))))
    )
    expect_false(
        any(is.nan(as.numeric(proxyC::simil(mat, method = "simple matching", margin = 2, use_nan = TRUE))))
    )

    # faith
    expect_false(
        any(is.nan(as.numeric(proxyC::simil(mat, method = "faith", margin = 1, use_nan = TRUE))))
    )
    expect_false(
        any(is.nan(as.numeric(proxyC::simil(mat, method = "faith", margin = 2, use_nan = TRUE))))
    )

})


test_that("simil works with dense matrices", {

    smat <- rsparsematrix(100, 50, 0.5)
    dmat <- as.matrix(smat)
    emat <- Matrix(smat, sparse = FALSE)
    s <- proxyC::simil(smat, smat)

    expect_identical(as.matrix(proxyC::simil(dmat, dmat)), as.matrix(s))
    expect_identical(as.matrix(proxyC::simil(emat, emat)), as.matrix(s))
    expect_silent(proxyC::simil(smat > 0, smat > 0))
    expect_error(proxyC::simil(forceSymmetric(emat), forceSymmetric(emat)))

})
