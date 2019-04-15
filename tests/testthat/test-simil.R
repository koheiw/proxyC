context("test simil")

test_mt <- Matrix::rsparsematrix(100, 100, 0.5)

test_simil <- function(x, method, margin, ignore_upper = FALSE, ignore_diag = TRUE, ...) {
    # test with only x
    s1 <- as.matrix(simil(x, method = method, margin = margin, ...))
    s2 <- proxy::as.matrix(proxy::simil(as.matrix(x),
                                method = method, by_rows = margin == 1, diag = TRUE, ...))
    if (ignore_diag)
        diag(s1) <- diag(s2) <- 0
    if (ignore_upper)
        s1[upper.tri(s1, TRUE)] <- s2[upper.tri(s2, TRUE)] <- 0
    expect_equal(as.numeric(s1), as.numeric(s2), tolerance = 0.001)

    # test with x and y, different size
    if (margin == 1) {
        y <- x[sample(nrow(x), 10),]
    } else {
        y <- x[,sample(ncol(x), 10)]
    }

    s3 <- as.matrix(simil(x, y, method = method, margin = margin, ...))
    s4 <- as.matrix(proxy::simil(as.matrix(x), as.matrix(y),
                                 method = method, by_rows = margin == 1, diag = TRUE, ...))
    if (ignore_diag)
        diag(s3) <- diag(s4) <- 0
    if (ignore_upper)
        s3[upper.tri(s3, TRUE)] <- s4[upper.tri(s4, TRUE)] <- 0
    expect_equal(as.numeric(s3), as.numeric(s4), tolerance = 0.001)


    # test with x and y, same size
    if (margin == 1) {
        y <- x[sample(nrow(x)),]
    } else {
        y <- x[,sample(ncol(x))]
    }

    s3 <- as.matrix(simil(x, y, method = method, margin = margin, ...))
    s4 <- as.matrix(proxy::simil(as.matrix(x), as.matrix(y),
                                 method = method, by_rows = margin == 1, diag = TRUE, ...))
    if (ignore_diag)
        diag(s3) <- diag(s4) <- 0
    if (ignore_upper)
        s3[upper.tri(s3, TRUE)] <- s4[upper.tri(s4, TRUE)] <- 0
    expect_equal(as.numeric(s3), as.numeric(s4), tolerance = 0.001)
}

test_that("test textstat_proxy cosine similarity", {
    skip_if_not_installed("proxy")
    test_simil(test_mt, "cosine", margin = 1)
    test_simil(test_mt, "cosine", margin = 2)
})

test_that("test textstat_proxy correlation similarity", {
    skip_if_not_installed("proxy")
    test_simil(test_mt, "correlation", margin = 1)
    test_simil(test_mt, "correlation", margin = 2)
})

test_that("test textstat_proxy jaccard similarity", {
    skip_if_not_installed("proxy")
    test_simil(test_mt, "jaccard", margin = 1)
    test_simil(test_mt, "jaccard", margin = 2)
})

test_that("test textstat_proxy ejaccard similarity", {
    skip_if_not_installed("proxy")
    test_simil(test_mt, "ejaccard", margin = 1)
    test_simil(test_mt, "ejaccard", margin = 2)
})

test_that("test textstat_proxy dice similarity", {
    skip_if_not_installed("proxy")
    test_simil(test_mt, "dice", margin = 1)
    test_simil(test_mt, "dice", margin = 2)
})

test_that("test textstat_proxy edice similarity", {
    skip_if_not_installed("proxy")
    test_simil(test_mt, "edice", margin = 1)
    test_simil(test_mt, "edice", margin = 2)
})

test_that("test textstat_proxy simple matching similarity", {
    skip_if_not_installed("proxy")
    test_simil(test_mt, "simple matching", margin = 1)
    test_simil(test_mt, "simple matching", margin = 2)
})

test_that("test textstat_proxy hamman similarity", {
    skip_if_not_installed("proxy")
    test_simil(test_mt, "hamman", margin = 1)
    test_simil(test_mt, "hamman", margin = 2)
})

test_that("test textstat_proxy faith similarity", {
    skip_if_not_installed("proxy")
    test_simil(test_mt, "faith", margin = 1, ignore_upper = TRUE)
    test_simil(test_mt, "faith", margin = 2, ignore_upper = TRUE)
})
