
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

    # ignore use_nan
    expect_all_false(
        is.nan(rowSds(as.matrix(mt), use_nan = TRUE))
    )
    expect_all_false(
        is.nan(colSds(as.matrix(mt), use_nan = TRUE))
    )

    # with NA
    mt2 <- Matrix(c(0.1, 0.2, NA, -0.1, -0.2, 0), nrow = 2, sparse = TRUE)
    expect_equal(
        is.nan(rowSds(mt2)),
        c(TRUE, FALSE)
    )
    expect_equal(
        is.nan(colSds(mt2)),
        c(FALSE, TRUE, FALSE)
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

    # ignore use_nan
    expect_all_false(
        is.nan(colMaxs(as.matrix(mt), use_nan = TRUE))
    )
    expect_all_false(
        is.nan(rowMaxs(as.matrix(mt), use_nan = TRUE))
    )
    expect_all_false(
        is.nan(colMins(as.matrix(mt), use_nan = TRUE))
    )
    expect_all_false(
        is.nan(rowMins(as.matrix(mt), use_nan = TRUE))
    )

    # with NA
    mt2 <- Matrix(c(0.1, 0.2, NA, -0.1, -0.2, 0), nrow = 2, sparse = TRUE)
    expect_equal(
        is.nan(colMaxs(mt2)),
        c(FALSE, FALSE, FALSE)
    )
    expect_equal(
        is.nan(rowMaxs(mt2)),
        c(FALSE, FALSE)
    )
    expect_equal(
        is.nan(colMins(mt2)),
        c(FALSE, FALSE, FALSE)
    )
    expect_equal(
        is.nan(rowMins(mt2)),
        c(FALSE, FALSE)
    )
})

test_that("maxRow and maxCol are working", {

    mt <- rsparsematrix(100, 100, 0.5)
    rownames(mt) <- paste0("row", 1:100)
    colnames(mt) <- paste0("col", 1:100)

    expect_identical(
        max.col(mt, ties.method = "first"),
        unname(maxCol(mt))
    )
    expect_identical(
        max.col(t(mt), ties.method = "first"),
        unname(maxRow(mt))
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

    # ignore use_nan
    expect_all_false(
        is.na(maxCol(as.matrix(mt), use_nan = TRUE))
    )
    expect_all_false(
        is.nan(maxRow(as.matrix(mt), use_nan = TRUE))
    )

    # with NA
    mt2 <- Matrix(c(0.1, 0.2, NA, -0.1, -0.2, 0), nrow = 2, sparse = TRUE)
    expect_equal(
        is.nan(maxCol(mt2)),
        c(FALSE, FALSE)
    )
    expect_equal(
        is.nan(maxRow(mt2)),
        c(FALSE, FALSE, FALSE)
    )
})
