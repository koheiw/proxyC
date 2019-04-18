context("test dist")

test_mt <- Matrix::rsparsematrix(100, 100, 0.5)

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

    s3 <- as.matrix(dist(x, y, method = method, margin = margin, ...))
    s4 <- as.matrix(proxy::dist(as.matrix(x), as.matrix(y),
                                method = method, by_rows = margin == 1, diag = TRUE, ...))
    if (ignore_upper)
        s3[upper.tri(s3, TRUE)] <- s4[upper.tri(s4, TRUE)] <- 0
    expect_equal(as.numeric(s3), as.numeric(s4), tolerance = 0.001)
}

test_that("test dist euclidean distance", {
    skip_if_not_installed("proxy")
    test_dist(test_mt, "euclidean", margin = 1)
    test_dist(test_mt, "euclidean", margin = 2)
})

# test_that("test kullback kullback distance", {
#     skip_if_not_installed("proxy")
#     # make dense matrix to avoide Inf in proxy::dist
#     test_mt_dense <- test_mt + 1
#     # proxy::dist() also incorrectly produces symmetric matrix
#     test_dist(test_mt_dense, "kullback", margin = 1, ignore_upper = TRUE)
#     test_dist(test_mt_dense, "kullback", margin = 2, ignore_upper = TRUE)
# })

test_that("test dist manhattan distance", {
    skip_if_not_installed("proxy")
    test_dist(test_mt, "manhattan", margin = 1)
    test_dist(test_mt, "manhattan", margin = 2)
})

test_that("test dist maximum distance", {
    skip_if_not_installed("proxy")
    test_dist(test_mt, "maximum", margin = 1)
    test_dist(test_mt, "maximum", margin = 2)
})

# test_that("test dist canberra distance", {
#     skip_if_not_installed("proxy")
#     test_dist(test_mt, "canberra", margin = 1)
#     test_dist(test_mt, "canberra", margin = 2)
# })

test_that("test dist minkowski distance", {
    skip_if_not_installed("proxy")
    test_dist(test_mt, "minkowski", margin = 1, p = 0.1)
    test_dist(test_mt, "minkowski", margin = 2, p = 0.1)
    test_dist(test_mt, "minkowski", margin = 1, p = 2)
    test_dist(test_mt, "minkowski", margin = 2, p = 2)
    test_dist(test_mt, "minkowski", margin = 1, p = 10)
    test_dist(test_mt, "minkowski", margin = 2, p = 10)
})
