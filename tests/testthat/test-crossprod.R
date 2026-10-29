
test_that("crossprod works", {

    mat1 <- Matrix::rsparsematrix(500, 100, 0.5)
    mat2 <- Matrix::rsparsematrix(500, 1000, 0.5)

    prod0 <- Matrix::crossprod(mat1, mat2)
    prod1 <- proxyC::crossprod(mat1, mat2)

    expect_identical(
        dim(prod1), dim(prod0)
    )

    expect_equal(
        as.matrix(prod1), as.matrix(prod0),
        tol = 1e-10
    )

    # with min_prod
    prod1_min <- proxyC::crossprod(mat1, mat2, min_prod = 1.0)
    prod0_min <- prod0
    prod0_min[prod0_min < 1.0] <- 0.0

    expect_equal(
        as.matrix(prod1_min), as.matrix(prod0_min),
        tol = 1e-10
    )

    prod2 <- proxyC::crossprod(mat1, mat2, sparse = FALSE)

    expect_true(
        is(prod2, "dgeMatrix")
    )
    expect_equal(
        as.matrix(prod2),
        as.matrix(prod0)
    )

    # large values

    mat10 <- mat1 * 10000
    mat20 <- mat2 * 10000

    prod00 <- Matrix::crossprod(mat10, mat20)
    prod10 <- proxyC::crossprod(mat10, mat20)

    expect_identical(
        dim(prod10), dim(prod00)
    )

    expect_equal(
        as.matrix(prod10), as.matrix(prod00),
        tol = 1e-10
    )

})


test_that("test tcrossprod", {

    mat1 <- Matrix::rsparsematrix(100, 500, 0.5)
    mat2 <- Matrix::rsparsematrix(1000, 500, 0.5)

    prod0 <- Matrix::tcrossprod(mat1, mat2)
    prod1 <- proxyC::tcrossprod(mat1, mat2)

    expect_identical(
        dim(prod1), dim(prod0)
    )

    expect_equal(
        as.matrix(prod1), as.matrix(prod0),
        tol = 1e-10
    )

    # with min_prod
    prod1_min <- proxyC::tcrossprod(mat1, mat2, min_prod = 1.0)
    prod0_min <- prod0
    prod0_min[prod0_min < 1.0] <- 0.0

    expect_equal(
        as.matrix(prod1_min), as.matrix(prod0_min),
        tol = 1e-10
    )

    prod2 <- proxyC::tcrossprod(mat1, mat2, sparse = FALSE)

    expect_true(
        is(prod2, "dgeMatrix")
    )
    expect_equal(
        as.matrix(prod2),
        as.matrix(prod0)
    )

})
