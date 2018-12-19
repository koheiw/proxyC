context("test proxy")

test_mt <- Matrix::rsparsematrix(100, 100, 0.5)

test_that("porxy stops as expected for methods not supported",{
    expect_error(simil(test_mt, method = "Yule"))
    expect_error(dist(test_mt, method = "Yule"))
    expect_error(proxyC:::proxy(test_mt, method = "Yule"))
})

test_that("raises error when p is smaller than 1", {
    expect_error(proxyC:::proxy(test_mt, method = "minkowski", p = 0))
    expect_error(proxyC:::proxy(test_mt, method = "minkowski", p = -1))
})

test_that("sparse objects are of expected class and occur when expected", {

    expect_is(proxyC:::proxy(test_mt),
              "dsTMatrix")
    expect_is(proxyC:::proxy(test_mt, min_proxy = 10),
              "dsTMatrix")
    expect_is(proxyC:::proxy(test_mt, rank = 2),
              "dgTMatrix")
    expect_is(proxyC:::proxy(test_mt, method = "kullback"),
              "dgTMatrix")

})

test_that("rank argument is working", {

    expect_error(proxyC:::proxy(test_mt, rank = 0),
                 "rank must be great than or equal to 1")

    expect_equal(as.matrix(proxyC:::proxy(test_mt)),
                 as.matrix(proxyC:::proxy(test_mt, rank = 100)))

    expect_equal(as.matrix(proxyC:::proxy(test_mt, rank = 3)),
                 apply(as.matrix(proxyC:::proxy(test_mt)), 2,
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

test_that("proxyC:::proxy raises error when the numnber of columns/rows are different", {
    expect_silent(
        proxyC:::proxy(test_mt[1:5,], test_mt[1:5,], margin = 2)
    )
    expect_error(proxyC:::proxy(test_mt[,1:5], test_mt[,1:10], margin = 1),
                 "x and y must have the same number of columns")
    expect_error(proxyC:::proxy(test_mt[1:5,], test_mt[1:10,], margin = 2),
                 "x and y must have the same number of rows")
})

test_that("proxyC:::proxy raises error when x or y is not a sparse matrix", {
    expect_error(proxyC:::proxy(as(test_mt, "denseMatrix")),
                 "x must be a sparseMatrix")
    expect_error(proxyC:::proxy(test_mt, as(test_mt, "denseMatrix")),
                 "y must be a sparseMatrix")
    expect_error(proxyC:::proxy(as.matrix(test_mt)),
                 "x must be a sparseMatrix")
    expect_error(proxyC:::proxy(test_mt, as.matrix(test_mt)),
                 "y must be a sparseMatrix")
})
