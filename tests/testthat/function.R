require(Matrix)

mat_test <- rsparsematrix(100, 50, 0.5)
mat_test[1, ] <- 0.0 # set sum(x) == 0
mat_test[2, ] <- 0.5 # set sd(x) == 0
mat_test[3, ] <- 0.0001 # set sd(x) == 0

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

    method_nan <- c("cosine", "correlation") # methods affected by all = 0 or sd = 0

    # test with only x
    suppressWarnings({
        s1 <- as.matrix(simil(x, method = method, margin = margin, use_nan = use_nan, ...))
    })
    s2 <- proxy::as.matrix(proxy::simil(as.matrix(x),
                                        method = method, by_rows = margin == 1, diag = TRUE, ...))

    if (use_nan) {
        if (method %in% method_nan) {
            s2[is_all0(x, x, margin = margin)] <- NaN
            if (method == "correlation")
                s2[is_sd0(x, x, margin = margin)] <- NaN
        }
    }
    if (ignore_diag)
        diag(s1) <- diag(s2) <- 0
    if (ignore_upper)
        s1[upper.tri(s1, TRUE)] <- s2[upper.tri(s2, TRUE)] <- 0

    expect_equal(as.numeric(s1), as.numeric(s2), tolerance = 0.001)

    # test with x and y, different size
    if (margin == 1) {
        y <- x[sample(nrow(x), pmin(nrow(x), 10)),]
    } else {
        y <- x[,sample(ncol(x), pmin(ncol(x), 10))]
    }

    suppressWarnings({
        s3 <- as.matrix(simil(x, y, method = method, margin = margin, use_nan = use_nan, ...))
    })
    s4 <- as.matrix(proxy::simil(as.matrix(x), as.matrix(y),
                                 method = method, by_rows = margin == 1, diag = TRUE, ...))

    if (use_nan) {
        if (method %in% method_nan) {
            s4[is_all0(x, y, margin = margin)] <- NaN
            if (method == "correlation")
                s4[is_sd0(x, y, margin = margin)] <- NaN
        }
    }
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

    suppressWarnings({
        s5 <- as.matrix(simil(x, y, method = method, margin = margin, use_nan = use_nan, ...))
    })
    s6 <- as.matrix(proxy::simil(as.matrix(x), as.matrix(y),
                                 method = method, by_rows = margin == 1, diag = TRUE, ...))
    if (use_nan) {
        if (method %in% method_nan) {
            s6[is_all0(x, y, margin = margin)] <- NaN
            if (method == "correlation")
                s6[is_sd0(x, y, margin = margin)] <- NaN
        }
    }
    if (ignore_diag)
        diag(s5) <- diag(s6) <- 0
    if (ignore_upper)
        s5[upper.tri(s5, TRUE)] <- s6[upper.tri(s6, TRUE)] <- 0
    expect_equal(as.numeric(s5), as.numeric(s6), tolerance = 0.001)
}


test_dist <- function(x, method, margin, ignore_upper = FALSE, ...) {

    # test with only x
    suppressWarnings({
        s1 <- as.matrix(dist(x, method = method, margin = margin, ...))
    })
    s2 <- as.matrix(proxy::dist(as.matrix(x),
                                method = method, by_rows = margin == 1, diag = TRUE, ...))

    if (ignore_upper)
        s1[upper.tri(s1, TRUE)] <- s2[upper.tri(s2, TRUE)] <- 0
    expect_equal(as.numeric(s1), as.numeric(s2), tolerance = 0.001)

    # test with x and y, different size
    if (margin == 1) {
        y <- x[sample(nrow(x), pmin(nrow(x), 10)),]
    } else {
        y <- x[,sample(ncol(x), pmin(ncol(x), 10))]
    }

    suppressWarnings({
        s3 <- as.matrix(dist(x, y, method = method, margin = margin, ...))
    })
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

    suppressWarnings({
        s5 <- as.matrix(dist(x, y, method = method, margin = margin, ...))
    })
    s6 <- as.matrix(proxy::dist(as.matrix(x), as.matrix(y),
                                method = method, by_rows = margin == 1, diag = TRUE, ...))
    if (ignore_upper)
        s5[upper.tri(s5, TRUE)] <- s6[upper.tri(s5, TRUE)] <- 0
    expect_equal(as.numeric(s5), as.numeric(s6), tolerance = 0.001)
}

