source("function.R")

test_that("test euclidean distance", {
    skip_if_not_installed("proxy")
    test_dist(mat_test, "euclidean", margin = 1)
    test_dist(mat_test, "euclidean", margin = 2)
})

test_that("test manhattan distance", {
    skip_if_not_installed("proxy")
    test_dist(mat_test, "manhattan", margin = 1)
    test_dist(mat_test, "manhattan", margin = 2)
})

test_that("test maximum distance", {
    skip_if_not_installed("proxy")
    test_dist(mat_test, "maximum", margin = 1)
    test_dist(mat_test, "maximum", margin = 2)
})

test_that("test minkowski distance", {
    skip_if_not_installed("proxy")
    test_dist(mat_test, "minkowski", margin = 1, p = 0.1)
    test_dist(mat_test, "minkowski", margin = 2, p = 0.1)
    test_dist(mat_test, "minkowski", margin = 1, p = 2)
    test_dist(mat_test, "minkowski", margin = 2, p = 2)
    test_dist(mat_test, "minkowski", margin = 1, p = 10)
    test_dist(mat_test, "minkowski", margin = 2, p = 10)
})

test_that("test canberra distance", {
    skip_if_not_installed("proxy")
    # proxyC and proxy disagree when sparsity is high
    smat <- rsparsematrix(100, 100, 0.99, rand.x = sample.int)
    test_dist(smat, "canberra", margin = 1)
    test_dist(smat, "canberra", margin = 2)
})

test_that("test chisquared distance", {
    skip_if_not_installed("entropy")

    smat <- rsparsematrix(10, 2, 0.5, rand.x = sample.int)
    expect_equal(
        proxyC::dist(smat, method = "chisquared", margin = 2)[1,2],
        0.0
    )
    dmat <- as.matrix(smat)
    expect_equal(
        proxyC::dist(smat, method = "chisquared", margin = 2, smooth = 1)[1,2],
        entropy::chi2indep.empirical(dmat[,c(1, 2)] + 1)
    )
})

test_that("test kullback leibler distance", {
    skip_if_not_installed("entropy")
    smat <- rsparsematrix(10, 2, 0.5, rand.x = sample.int)
    expect_equal(
        proxyC::dist(smat, method = "kullback", margin = 2)[1,2],
        0.0
    )
    dmat <- as.matrix(smat)
    expect_equal(
        as.matrix(proxyC::dist(smat, method = "kullback", margin = 2, smooth = 1))[1,2],
        entropy::KL.empirical(dmat[,1] + 1, dmat[,2] + 1)
    )
    expect_equal(
        as.matrix(proxyC::dist(smat, method = "kullback", margin = 2, smooth = 1))[2,1],
        entropy::KL.empirical(dmat[,2] + 1, dmat[,1] + 1)
    )
})

test_that("test jeffreys distance", {
    skip_if_not_installed("entropy")
    smat <- rsparsematrix(10, 2, 0.5, rand.x = sample.int)
    expect_equal(
        proxyC::dist(smat, method = "jeffreys", margin = 2)[1,2],
        0.0
    )
    dmat <- as.matrix(smat)
    kl <- as.matrix(proxyC::dist(smat, method = "kullback", margin = 2, smooth = 1))
    jd <- as.matrix(proxyC::dist(smat, method = "jeffreys", margin = 2, smooth = 1))
    dimnames(kl) <- dimnames(jd) <- list(NULL, NULL)
    expect_equal(kl + t(kl), jd)
})

test_that("test hamming distance", {
    new_mat_test <- rsparsematrix(100, 90, 1, rand.x = function(x) sample.int(10, x, replace = TRUE))
    dmat <- as.matrix(proxyC::dist(new_mat_test, method = "hamming"))
    dmat_manual <-
        sapply(seq_len(nrow(new_mat_test)), function(i) {
            rowSums(sweep(new_mat_test, 2, new_mat_test[i, ], "!="))
        })
    expect_equal(
        dmat,
        dmat_manual,
        check.attributes = FALSE
    )
    expect_equal(
        mean(dmat[!diag(nrow(dmat))]),
        .9 * nrow(new_mat_test), # thanks to rand.x function, there's a 10% chance that values from different rows will match
        tolerance = 1
    )
})

test_that("use_nan is working", {

    mat1 <- Matrix::Matrix(1:4, nrow = 1, sparse = TRUE)
    mat2 <- Matrix::Matrix(rep(0, 4), nrow = 1, sparse = TRUE)

    expect_warning(proxyC::dist(mat1, mat2, method = "kullback", use_nan = FALSE),
                   "x or y has vectors with all zero; consider setting use_nan = TRUE")
    expect_warning(proxyC::dist(mat1, mat2, method = "chisquared", use_nan = FALSE),
                   "x or y has vectors with all zero; consider setting use_nan = TRUE")

    suppressWarnings({
        expect_equal(proxyC::dist(mat1, mat2, method = "kullback", use_nan = FALSE)[1,1], 0)
        expect_equal(proxyC::dist(mat1, mat2, method = "kullback", use_nan = FALSE, diag = TRUE)[1,1], 0)
    })

    expect_true(is.na(proxyC::dist(mat1, mat2, method = "kullback", use_nan = TRUE)[1,1]))
    expect_true(is.na(proxyC::dist(mat1, mat2, method = "kullback", use_nan = TRUE, diag = TRUE)[1,1]))

    suppressWarnings({
        expect_equal(proxyC::dist(mat1, mat2, method = "chisquared", use_nan = FALSE)[1,1], 0)
        expect_equal(proxyC::dist(mat1, mat2, method = "chisquared", use_nan = FALSE, diag = TRUE)[1,1], 0)
    })

    expect_true(is.na(proxyC::dist(mat1, mat2, method = "chisquared", use_nan = TRUE)[1,1]))
    expect_true(is.na(proxyC::dist(mat1, mat2, method = "chisquared", use_nan = TRUE, diag = TRUE)[1,1]))
})

test_that("dist returns zero or NaN correctly", {

    mat <- Matrix::Matrix(matrix(c(0, 0, 0,
                                   1, 1, 1,
                                   1, 5, 2,
                                   2, 3, 4), byrow = TRUE, nrow = 4), sparse = TRUE)

    # euclidean
    expect_equivalent(
        as.matrix(proxyC::dist(mat, method = "euclidean", margin = 1, use_nan = TRUE) == 0),
        as.matrix(bandSparse(4, 4, k = 0))
    )
    expect_equivalent(
        as.matrix(proxyC::dist(mat, method = "euclidean", margin = 2, use_nan = TRUE) == 0),
        as.matrix(bandSparse(3, 3, k = 0))
    )

    # kullback
    expect_equivalent(
        suppressWarnings(as.matrix(proxyC::dist(mat, method = "kullback", margin = 1, use_nan = FALSE) == 0)),
        is_all0(mat, margin = 1) | as.matrix(bandSparse(4, 4, k = 0))
    )
    expect_equivalent(
        suppressWarnings(as.matrix(proxyC::dist(mat, method = "kullback", margin = 2, use_nan = FALSE) == 0)),
        matrix(TRUE, 3, 3)
    )
    expect_equivalent(
        is.nan(as.matrix(proxyC::dist(mat, method = "kullback", margin = 1, use_nan = TRUE))),
        is_all0(mat, margin = 1)
    )
    expect_equivalent(
        is.nan(as.matrix(proxyC::dist(mat, method = "kullback", margin = 2, use_nan = TRUE))),
        matrix(TRUE, 3, 3)
    )
    expect_equivalent(
        is.nan(as.matrix(proxyC::dist(mat, method = "kullback", margin = 1, smooth = 1, use_nan = TRUE))),
        matrix(FALSE, 4, 4)
    )
    expect_equivalent(
        is.nan(as.matrix(proxyC::dist(mat, method = "kullback", margin = 2, smooth = 1, use_nan = TRUE))),
        matrix(FALSE, 3, 3)
    )

    # chisquared
    expect_equivalent(
        suppressWarnings(as.matrix(proxyC::dist(mat, method = "chisquared", margin = 1, use_nan = FALSE) == 0)),
        is_all0(mat, margin = 1) | as.matrix(bandSparse(4, 4, k = 0))
    )
    expect_equivalent(
        suppressWarnings(as.matrix(proxyC::dist(mat, method = "chisquared", margin = 2, use_nan = FALSE) == 0)),
        matrix(TRUE, 3, 3)
    )
    expect_equivalent(
        is.nan(as.matrix(proxyC::dist(mat, method = "chisquared", margin = 1, use_nan = TRUE))),
        is_all0(mat, margin = 1)
    )
    expect_equivalent(
        is.nan(as.matrix(proxyC::dist(mat, method = "chisquared", margin = 2, use_nan = TRUE))),
        matrix(TRUE, 3, 3)
    )
    expect_equivalent(
        is.nan(as.matrix(proxyC::dist(mat, method = "chisquared", margin = 1, smooth = 1, use_nan = TRUE))),
        matrix(FALSE, 4, 4)
    )
    expect_equivalent(
        is.nan(as.matrix(proxyC::dist(mat, method = "chisquared", margin = 2, smooth = 1, use_nan = TRUE))),
        matrix(FALSE, 3, 3)
    )

    # manhattan
    expect_false(
        any(is.nan(as.numeric(proxyC::dist(mat, method = "manhattan", margin = 1, use_nan = TRUE))))
    )
    expect_false(
        any(is.nan(as.numeric(proxyC::dist(mat, method = "manhattan", margin = 2, use_nan = TRUE))))
    )

    # maximum
    expect_false(
        any(is.nan(as.numeric(proxyC::dist(mat, method = "maximum", margin = 1, use_nan = TRUE))))
    )
    expect_false(
        any(is.nan(as.numeric(proxyC::dist(mat, method = "maximum", margin = 2, use_nan = TRUE))))
    )

    # canberra
    expect_false(
        any(is.nan(as.numeric(proxyC::dist(mat, method = "canberra", margin = 1, use_nan = TRUE))))
    )
    expect_false(
        any(is.nan(as.numeric(proxyC::dist(mat, method = "canberra", margin = 2, use_nan = TRUE))))
    )

    # minkowski
    expect_false(
        any(is.nan(as.numeric(proxyC::dist(mat, method = "minkowski", margin = 1, use_nan = TRUE))))
    )
    expect_false(
        any(is.nan(as.numeric(proxyC::dist(mat, method = "minkowski", margin = 2, use_nan = TRUE))))
    )

    # hamming
    expect_false(
        any(is.nan(as.numeric(proxyC::dist(mat, method = "hamming", margin = 1, use_nan = TRUE))))
    )
    expect_false(
        any(is.nan(as.numeric(proxyC::dist(mat, method = "hamming", margin = 2, use_nan = TRUE))))
    )

})


test_that("dist works with dense matrices", {

    smat <- rsparsematrix(100, 50, 0.5)
    dmat <- as.matrix(smat)
    emat <- Matrix(smat, sparse = FALSE)
    d <- proxyC::dist(smat, smat)

    expect_identical(as.matrix(proxyC::dist(dmat, dmat)), as.matrix(d))
    expect_identical(as.matrix(proxyC::dist(emat, emat)), as.matrix(d))
    expect_silent(proxyC::dist(smat > 0, smat > 0))
    expect_error(proxyC::dist(forceSymmetric(emat), forceSymmetric(emat)))

})
