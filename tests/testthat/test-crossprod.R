
test_that("test crossprod", {


    mat1 <- Matrix::rsparsematrix(500, 100, 0.5)
    mat2 <- Matrix::rsparsematrix(500, 1000, 0.5)

    prod1 <- proxyC:::crossprod(mat1, mat2)
    prod2 <- Matrix::crossprod(mat1, mat2)

    expect_identical(
        dim(prod1), dim(prod2)
    )

    expect_equal(
        as.matrix(prod1), as.matrix(prod2),
        tol = 1e-10
    )

    # with min_prod
    prod1_min <- proxyC:::crossprod(mat1, mat2, min_prod = 1.0)
    prod2_min <- prod2
    prod2_min[prod2_min <= 1.0] <- 0.0

    expect_equal(
        as.matrix(prod1_min), as.matrix(prod2_min),
        tol = 1e-10
    )

})
