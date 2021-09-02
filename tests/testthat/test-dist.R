require(Matrix)
mat_test <- rsparsematrix(100, 90, 0.5)
mat_test[1, ] <- 10.123 # add one row with sd(x) == 0
mat_test[, 1] <- 10.123 # add one col with sd(x) == 0

test_dist <- function(x, method, margin, ignore_upper = FALSE, ...) {
    # test with only x
    s1 <- as.matrix(dist(x, method = method, margin = margin, ...))
    s2 <- as.matrix(proxy::dist(as.matrix(x),
                                method = method, by_rows = margin == 1, diag = TRUE, ...))

    if (ignore_upper)
        s1[upper.tri(s1, TRUE)] <- s2[upper.tri(s2, TRUE)] <- 0
    expect_equal(as.numeric(s1), as.numeric(s2), tolerance = 0.001)

    # test with x and y, different size
    if (margin == 1) {
        y <- x[sample(nrow(x), 10),]
    } else {
        y <- x[,sample(ncol(x), 10)]
    }

    s3 <- as.matrix(dist(x, y, method = method, margin = margin, ...))
    s4 <- as.matrix(proxy::dist(as.matrix(x), as.matrix(y),
                                method = method, by_rows = margin == 1, diag = TRUE, ...))
    if (ignore_upper)
        s3[upper.tri(s3, TRUE)] <- s4[upper.tri(s4, TRUE)] <- 0
    expect_equal(as.numeric(s3), as.numeric(s4), tolerance = 0.001)

    # test with x and y, same size
    if (margin == 1) {
        y <- x[sample(nrow(x)),]
    } else {
        y <- x[,sample(ncol(x))]
    }

    s5 <- as.matrix(dist(x, y, method = method, margin = margin, ...))
    s6 <- as.matrix(proxy::dist(as.matrix(x), as.matrix(y),
                                method = method, by_rows = margin == 1, diag = TRUE, ...))
    if (ignore_upper)
        s5[upper.tri(s5, TRUE)] <- s6[upper.tri(s5, TRUE)] <- 0
    expect_equal(as.numeric(s5), as.numeric(s6), tolerance = 0.001)
}

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
