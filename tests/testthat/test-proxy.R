context("test proxy")

dfmt_test <- Matrix::rsparsematrix(100, 100, 0.5)

test_that("porxy takes only sparse matrix",{
    expect_error(simil(matrix(dfmt_test)), "x must be a sparseMatrix")
    expect_error(dist(matrix(dfmt_test)), "x must be a sparseMatrix")
})

test_that("porxy stops as expected for methods not supported",{
    expect_error(simil(dfmt_test, method = "Yule"))
    expect_error(dist(dfmt_test, method = "Yule"))
    expect_error(proxyC:::proxy(dfmt_test, method = "Yule"))
})

test_that("raises error when p is smaller than 1", {
    expect_error(proxyC:::proxy(dfmt_test, method = "minkowski", p = 0))
    expect_error(proxyC:::proxy(dfmt_test, method = "minkowski", p = -1))
})

test_that("sparse objects are of expected class and occur when expected", {

    expect_is(proxyC:::proxy(dfmt_test),
              "dsTMatrix")
    expect_is(proxyC:::proxy(dfmt_test, min_proxy = 10),
              "dsTMatrix")
    expect_is(proxyC:::proxy(dfmt_test, rank = 2),
              "dgTMatrix")
    expect_is(proxyC:::proxy(dfmt_test, method = "kullback"),
              "dgTMatrix")

})

test_that("rank argument is working", {

    expect_error(proxyC:::proxy(dfmt_test, rank = 0),
                 "rank must be great than or equal to 1")

    expect_equal(as.matrix(proxyC:::proxy(dfmt_test)),
                 as.matrix(proxyC:::proxy(dfmt_test, rank = 100)))

    expect_equal(as.matrix(proxyC:::proxy(dfmt_test, rank = 3)),
                 apply(as.matrix(proxyC:::proxy(dfmt_test)), 2,
                       function(x) ifelse(x >= sort(x, decreasing = TRUE)[3], x, 0)))

})

test_that("record zeros even in the sparse matrix", {
    mt <- Matrix::rsparsematrix(100, 100, 0.01)
    expect_true(any(proxyC:::proxy(mt)@x == 0))
    expect_true(any(proxyC:::proxy(mt, method = "cosine")@x == 0))
    expect_true(any(proxyC:::proxy(mt, method = "cosine", min_proxy = -0.5)@x == 0))
    expect_true(any(proxyC:::proxy(mt, method = "cosine", rank = 2)@x == 0))
    expect_true(any(proxyC:::proxy(mt, method = "dice")@x == 0))
})

test_that("do not record zeros when drop0 is TRUE", {
    mt <- Matrix::rsparsematrix(100, 100, 0.01)
    expect_false(any(proxyC:::proxy(mt, drop0 = TRUE)@x == 0))
    expect_false(any(proxyC:::proxy(mt, method = "cosine", drop0 = TRUE)@x == 0))
    expect_false(any(proxyC:::proxy(mt, method = "cosine", drop0 = TRUE, min_proxy = -0.5)@x == 0))
    expect_false(any(proxyC:::proxy(mt, method = "cosine", drop0 = TRUE, rank = 2)@x == 0))
    expect_false(any(proxyC:::proxy(mt, method = "dice", drop0 = TRUE)@x == 0))
})


test_that("proxyC:::proxy raises error when the numnber of columns/rows are different", {
    expect_silent(
        proxyC:::proxy(dfmt_test[1:5,], dfmt_test[1:5,], margin = 2)
    )
    expect_error(proxyC:::proxy(dfmt_test[,1:5], dfmt_test[,1:10], margin = 1),
                 "x and y must have the same number of columns")
    expect_error(proxyC:::proxy(dfmt_test[1:5,], dfmt_test[1:10,], margin = 2),
                 "x and y must have the same number of rows")
})

test_that("proxyC:::proxy raises error when x or y is not a sparse matrix", {
    expect_error(proxyC:::proxy(as(dfmt_test, "denseMatrix")),
                 "x must be a sparseMatrix")
    expect_error(proxyC:::proxy(dfmt_test, as(dfmt_test, "denseMatrix")),
                 "y must be a sparseMatrix")
    expect_error(proxyC:::proxy(as.matrix(dfmt_test)),
                 "x must be a sparseMatrix")
    expect_error(proxyC:::proxy(dfmt_test, as.matrix(dfmt_test)),
                 "y must be a sparseMatrix")
})

test_that("digits is working", {
    mat <- Matrix::Matrix(c(1, 2, 5, 3), ncol = 2, sparse = TRUE)
    sim1 <- simil(mat, method = "cosine")
    expect_true(all(sim1 <= 1))
    sim2 <- simil(mat, method = "cosine", digits = 3)
    expect_identical(as.numeric(as.matrix(sim2)), c(1.0, 0.925, 0.925, 1.0))
    sim3 <- simil(mat, method = "cosine", digits = 1)
    expect_identical(as.numeric(as.matrix(sim3)), c(1.0, 0.9, 0.9, 1.0))
})

test_that("colSds and rowSds are working", {
    mt <- Matrix::rsparsematrix(100, 100, 0.01)
    expect_equal(rowSds(mt), apply(mt, 1, sd))
    expect_equal(colSds(mt), apply(mt, 2, sd))
    expect_error(rowSds(matrix(mt)), "x must be a sparseMatrix")
    expect_error(colSds(matrix(mt)), "x must be a sparseMatrix")
})

test_that("colZeros and rowZeros are working", {
    mt <- Matrix::rsparsematrix(100, 100, 0.01)
    expect_equal(rowZeros(mt), apply(mt, 1, function(x) sum(x == 0)))
    expect_equal(colZeros(mt), apply(mt, 2, function(x) sum(x == 0)))
    expect_error(rowZeros(matrix(mt)), "x must be a sparseMatrix")
    expect_error(colZeros(matrix(mt)), "x must be a sparseMatrix")
})
