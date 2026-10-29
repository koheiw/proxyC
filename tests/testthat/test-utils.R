
test_that("colZeros and rowZeros are working", {

    mt <- rsparsematrix(100, 100, 0.01)
    dimnames(mt) <- list(paste0("row", seq_len(nrow(mt))),
                         paste0("col", seq_len(ncol(mt))))

    expect_equal(names(rowZeros(mt)), rownames(mt))
    expect_equal(names(colZeros(mt)), colnames(mt))
    expect_equal(names(rowSds(mt)), rownames(mt))
    expect_equal(names(colSds(mt)), colnames(mt))
    expect_equal(rowZeros(mt), apply(mt, 1, function(x) sum(x == 0)))
    expect_equal(colZeros(mt), apply(mt, 2, function(x) sum(x == 0)))
    expect_equal(rowZeros(as.matrix(mt)), apply(mt, 1, function(x) sum(x == 0)))
    expect_equal(colZeros(as.matrix(mt)), apply(mt, 2, function(x) sum(x == 0)))
})

test_that("colSds and rowSds are working", {

    mt <- rsparsematrix(100, 100, 0.01)
    dimnames(mt) <- list(paste0("row", seq_len(nrow(mt))),
                         paste0("col", seq_len(ncol(mt))))

    expect_equal(rowSds(mt), apply(mt, 1, sd))
    expect_equal(colSds(mt), apply(mt, 2, sd))
    expect_equal(rowSds(as.matrix(mt)), apply(mt, 1, sd))
    expect_equal(colSds(as.matrix(mt)), apply(mt, 2, sd))

    # enable use_nan
    expect_identical(
        is.nan(rowSds(mt, use_nan = TRUE)),
        rowZeros(mt) == 100
    )
    expect_identical(
        is.nan(colSds(mt, use_nan = TRUE)),
        colZeros(mt) == 100
    )

    # error
    expect_error(
        maxCol(rowSds(mt), use_nan = TRUE),
        "use_nan applies only x is a sparse matrix"
    )
    expect_error(
        maxRow(colSds(mt), use_nan = TRUE),
        "use_nan applies only x is a sparse matrix"
    )
})

test_that("colMax, rowMax, colMin, rowMin are working", {

    mt <- rsparsematrix(100, 100, 0.01)
    dimnames(mt) <- list(paste0("row", seq_len(nrow(mt))),
                         paste0("col", seq_len(ncol(mt))))


    expect_equal(colMaxs(mt), apply(mt, 2, max))
    expect_equal(names(colMaxs(mt)), colnames(mt))

    expect_equal(rowMaxs(mt), apply(mt, 1, max))
    expect_equal(names(rowMaxs(mt)), rownames(mt))

    expect_equal(colMins(mt), apply(mt, 2, min))
    expect_equal(names(colMaxs(mt)), colnames(mt))

    expect_equal(rowMins(mt), apply(mt, 1, min))
    expect_equal(names(rowMins(mt)), rownames(mt))

    # disable use_nan
    expect_all_false(
        is.nan(colMaxs(mt, use_nan = FALSE))
    )
    expect_all_false(
        is.nan(rowMaxs(mt, use_nan = FALSE))
    )
    expect_all_false(
        is.nan(colMins(mt, use_nan = FALSE))
    )
    expect_all_false(
        is.nan(rowMins(mt, use_nan = FALSE))
    )

    # enable use_nan
    expect_identical(
        is.nan(colMaxs(mt, use_nan = TRUE)),
        colZeros(mt) == 100
    )
    expect_identical(
        is.nan(rowMaxs(mt, use_nan = TRUE)),
        rowZeros(mt) == 100
    )
    expect_identical(
        is.nan(colMins(mt, use_nan = TRUE)),
        colZeros(mt) == 100
    )
    expect_identical(
        is.nan(rowMins(mt, use_nan = TRUE)),
        rowZeros(mt) == 100
    )

    # error
    expect_error(
        maxCol(colMaxs(mt), use_nan = TRUE),
        "use_nan applies only x is a sparse matrix"
    )
    expect_error(
        maxRow(rowMaxs(mt), use_nan = TRUE),
        "use_nan applies only x is a sparse matrix"
    )
    expect_error(
        maxCol(colMins(mt), use_nan = TRUE),
        "use_nan applies only x is a sparse matrix"
    )
    expect_error(
        maxRow(rowMins(mt), use_nan = TRUE),
        "use_nan applies only x is a sparse matrix"
    )
})

test_that("maxRow and maxCol works with sparse marix", {

    mt <- rsparsematrix(100, 100, 0.01)
    rownames(mt) <- paste0("row", 1:100)
    colnames(mt) <- paste0("col", 1:100)

    # random
    expect_false(
        all(maxCol(mt) == maxCol(mt, ties_method = "first"))
    )
    expect_false(
        all(maxCol(mt) == maxCol(mt, ties_method = "last"))
    )
    expect_false(
        all(maxCol(mt) == maxCol(mt))
    )
    expect_identical(
        {set.seed(1234); maxCol(mt)},
        {set.seed(1234); maxCol(mt)}
    )

    # first
    expect_identical(
        max.col(mt, ties.method = "first"),
        unname(maxCol(mt, ties_method = "first"))
    )
    expect_identical(
        max.col(t(mt), ties.method = "first"),
        unname(maxRow(mt, ties_method = "first"))
    )

    # last
    expect_identical(
        max.col(mt, ties.method = "last"),
        unname(maxCol(mt, ties_method = "last"))
    )
    expect_identical(
        max.col(t(mt), ties.method = "last"),
        unname(maxRow(mt, ties_method = "last"))
    )

    expect_identical(
        names(maxRow(mt)),
        colnames(mt)
    )
    expect_identical(
        names(maxCol(mt)),
        rownames(mt)
    )

    # enable use_nan
    mt2 <- rsparsematrix(100, 100, 0.01)
    expect_identical(
        is.na(maxCol(mt2, use_nan = TRUE)), # NOTE: change to use_na?
        rowZeros(mt2) == 100
    )
    expect_identical(
        is.na(maxRow(mt2, use_nan = TRUE)),
        colZeros(mt2) == 100
    )

    # error
    expect_error(
        maxRow(mt, use_nan = c(TRUE, FALSE)),
        "use_nan must a single logical value"
    )

    # extreme case
    mt2 <- Matrix::Matrix(rep(NA_real_, 6), nrow = 2, sparse = TRUE)
    expect_equal(
        maxCol(mt2),
        rep(NA_integer_, 2)
    )
    expect_equal(
        maxRow(mt2),
        rep(NA_integer_, 3)
    )
})

test_that("maxRow and maxCol works with dense marix", {

    mt <- as.matrix(rsparsematrix(100, 100, 0.01))
    rownames(mt) <- paste0("row", 1:100)
    colnames(mt) <- paste0("col", 1:100)

    # random
    expect_false(
        all(maxCol(mt) == maxCol(mt, ties_method = "first"))
    )
    expect_false(
        all(maxCol(mt) == maxCol(mt, ties_method = "last"))
    )
    expect_false(
        all(maxCol(mt) == maxCol(mt))
    )
    expect_identical(
        {set.seed(1234); maxCol(mt)},
        {set.seed(1234); maxCol(mt)}
    )

    # first
    expect_identical(
        max.col(mt, ties.method = "first"),
        unname(maxCol(mt, ties_method = "first"))
    )
    expect_identical(
        max.col(t(mt), ties.method = "first"),
        unname(maxRow(mt, ties_method = "first"))
    )

    # last
    expect_identical(
        max.col(mt, ties.method = "last"),
        unname(maxCol(mt, ties_method = "last"))
    )
    expect_identical(
        max.col(t(mt), ties.method = "last"),
        unname(maxRow(mt, ties_method = "last"))
    )

    expect_identical(
        names(maxRow(mt)),
        colnames(mt)
    )
    expect_identical(
        names(maxCol(mt)),
        rownames(mt)
    )

    # error
    expect_error(
        maxCol(mt, use_nan = TRUE),
        "use_nan applies only x is a sparse matrix"
    )
    expect_error(
        maxRow(mt, use_nan = TRUE),
        "use_nan applies only x is a sparse matrix"
    )

    # extreme case
    mt2 <- matrix(rep(NA_real_, 6), nrow = 2)
    expect_equal(
        maxCol(mt2),
        rep(NA_integer_, 2)
    )
    expect_equal(
        maxRow(mt2),
        rep(NA_integer_, 3)
    )
})

test_that("handle NA properly", {

    mt <- Matrix(c(0.1, 0.2, NA, -0.1, -0.2, 0),
                  nrow = 2, sparse = TRUE)

    # colSds/rowSds
    expect_equal(
        is.nan(colSds(mt)),
        c(FALSE, TRUE, FALSE)
    )
    expect_equal(
        is.nan(rowSds(mt)),
        c(TRUE, FALSE)
    )

    # colMaxs/rowMaxs
    expect_equal(
        is.nan(colMaxs(mt)),
        c(FALSE, FALSE, FALSE)
    )
    expect_equal(
        is.nan(rowMaxs(mt)),
        c(FALSE, FALSE)
    )
    expect_equal(
        is.nan(colMins(mt)),
        c(FALSE, FALSE, FALSE)
    )
    expect_equal(
        is.nan(rowMins(mt)),
        c(FALSE, FALSE)
    )

    # maxCol/maxRow
    expect_equal(
        is.nan(maxCol(mt)),
        c(FALSE, FALSE)
    )
    expect_equal(
        is.nan(maxRow(mt)),
        c(FALSE, FALSE, FALSE)
    )
})

test_that("convertMatrix works", {

    mat100 <- rsparsematrix(100, 100, 1.0)
    mat50 <- rsparsematrix(100, 100, 0.5)
    mat1 <- rsparsematrix(100, 100, 0.01)
    mat0 <- rsparsematrix(100, 100, 0)

    expect_s4_class(
        proxyC:::convertMatrix(mat100),
        "dgCMatrix"
    )
    expect_s4_class(
        proxyC:::convertMatrix(mat50),
        "dgCMatrix"
    )
    expect_s4_class(
        proxyC:::convertMatrix(mat1),
        "dgCMatrix"
    )
    expect_s4_class(
        proxyC:::convertMatrix(mat0),
        "dgCMatrix"
    )

    expect_is(
        proxyC:::convertMatrix(as.matrix(mat100)),
        "matrix"
    )
    expect_is(
        proxyC:::convertMatrix(as.matrix(mat100)),
        "matrix"
    )
    expect_is(
        proxyC:::convertMatrix(as.matrix(mat1)),
        "matrix"
    )
    expect_is(
        proxyC:::convertMatrix(as.matrix(mat0)),
        "matrix"
    )
})

test_that("functions works with empty matrices", {

    smat <- rsparsematrix(0, 10, 0)
    dmat <- as.matrix(smat)

    expect_identical(rowSds(smat), numeric())
    expect_identical(colSds(smat), rep(NaN, 10))
    expect_identical(rowSds(dmat), numeric())
    expect_identical(colSds(dmat), rep(NaN, 10))

    expect_identical(rowZeros(smat), integer())
    expect_identical(colZeros(smat), integer(10))
    expect_identical(rowZeros(dmat), integer())
    expect_identical(colZeros(dmat), integer(10))

    expect_identical(rowMins(smat), numeric())
    expect_identical(colMins(smat), rep(Inf, 10))
    expect_identical(rowMins(dmat), numeric())
    expect_identical(colMins(dmat), rep(Inf, 10))

    expect_identical(rowMaxs(smat), numeric())
    expect_identical(colMaxs(smat), rep(-Inf, 10))
    expect_identical(rowMaxs(dmat), numeric())
    expect_identical(colMaxs(dmat), rep(-Inf, 10))

    expect_identical(maxRow(smat), rep(NA_integer_, 10))
    expect_identical(maxCol(smat), integer())
    expect_identical(maxRow(dmat), rep(NA_integer_, 10))
    expect_identical(maxCol(dmat), integer())

})
