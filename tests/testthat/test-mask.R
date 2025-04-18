require(Matrix)

set.seed(1234)
mat1_test <- rsparsematrix(50, 20, 0.5)
colnames(mat1_test) <- rep(c("a", "b", "c", "d", "e"), each = 4)
mat2_test <- rsparsematrix(50, 5, 0.5)
colnames(mat2_test) <- c("a", "b", "c", "d", "e")
msk_test <- mask(colnames(mat1_test), colnames(mat2_test))

test_that("mask is working", {

    # characeter vectors
    c1 <- rep(c("a", "b", "c", "d", "e"), each = 4)
    c2 <- c("a", "b", "c", "d", "e")

    msk1 <- mask(c1, c2)
    expect_equivalent(class(msk1), "lgTMatrix")
    expect_identical(rownames(msk1), c1)
    expect_identical(colnames(msk1), c2)

    expect_equal(
        dim(mask(character(), character())),
        c(0, 0)
    )

    expect_true(
        isSymmetric(mask(c1))
    )

    expect_error(
        mask(c1, as.factor(c2)),
        "x and y must be the same type of vectors"
    )

    # numeric vectors
    n1 <- rep(1:5, each = 4)
    n2 <- 1:5

    msk2 <- mask(n1, n2)
    expect_equivalent(class(msk2), "lgTMatrix")
    expect_null(rownames(msk2))
    expect_null(colnames(msk2))

    expect_equal(
        dim(mask(numeric(), numeric())),
        c(0, 0)
    )

    expect_true(
        isSymmetric(mask(n1))
    )

    expect_error(
        mask(n1, as.character(n2)),
        "x and y must be the same type of vectors"
    )
})

test_that("mask works with simil linear", {

    # masked dense
    s1 <- simil(mat1_test, mat2_test, margin = 2,
                mask = msk_test, use_nan = TRUE)
    expect_identical(as.matrix(msk_test != 0), as.matrix(!is.na(s1)))
    expect_identical(as.matrix(msk_test == 0), as.matrix(is.na(s1)))

    # masked ranked
    s2 <- simil(mat1_test, mat2_test, margin = 2, rank = 1, use_nan = FALSE,
                mask = msk_test)
    expect_true(all(colSums(s2 != 0) == 1))
    expect_identical(colSums(s2), sapply(split(s2@x, rownames(s2)[s2@i + 1]), max))

    # masked min
    s3 <- simil(mat1_test, mat2_test, margin = 2, min_simil = 0, use_nan = FALSE,
                mask = msk_test)
    expect_true(all(colSums(s3 != 0) >= 1))
    expect_true(all(sapply(split(s3@x, rownames(s3)[s3@i + 1]), min) > 0))

    # masked min with nan
    s4 <- simil(mat1_test, mat2_test, margin = 2, min_simil = -0.1, use_nan = TRUE,
                mask = msk_test)
    expect_true(all(colSums(apply(s4, 2, is.na)) >= 1))
    expect_true(all(sapply(split(s4@x, rownames(s4)[s4@i + 1]), min, na.rm = TRUE) > -0.1))

    # symmetric
    s5 <- simil(mat1_test, margin = 2, use_nan = FALSE,
                mask = mask(colnames(mat1_test)))
    expect_true(isSymmetric(s5))

    expect_error(
        simil(mat1_test, mat2_test, margin = 2, mask = msk_test[,-1]),
        "The shape of mask must be 20 x 5"
    )
})

test_that("mask works with simil pair", {

    # masked dense
    s1 <- simil(mat1_test, mat2_test, margin = 2, method = "jaccard",
                mask = msk_test)
    expect_identical(as.matrix(msk_test != 0), as.matrix(s1 != 0))

    # masked ranked
    s2 <- simil(mat1_test, mat2_test, margin = 2, method = "jaccard", rank = 1, drop0 = TRUE,
                mask = msk_test)
    expect_identical(colSums(s2), sapply(split(s2@x, rownames(s2)[s2@i + 1]), max))

    # masked min
    s3 <- simil(mat1_test, mat2_test, margin = 2,method = "jaccard",  min_simil = 0, drop0 = TRUE,
                mask = msk_test)
    expect_true(all(sapply(split(s3@x, rownames(s3)[s3@i + 1]), min) > 0))

    # masked min
    s4 <- simil(mat1_test, mat2_test, margin = 2, method = "jaccard", min_simil = -0.1, drop0 = TRUE,
                mask = msk_test)
    expect_true(all(sapply(split(s4@x, rownames(s4)[s4@i + 1]), min) > -0.1))

    expect_error(
        simil(mat1_test, mat2_test, margin = 2, method = "jaccard", mask = msk_test[,-1]),
        "The shape of mask must be 20 x 5"
    )
})

test_that("mask = NULL is the same as all TRUE", {

    s1 <- simil(mat1_test, mat2_test, margin = 2,
                mask = NULL)

    s2 <- simil(mat1_test, mat2_test, margin = 2,
                mask = Matrix(TRUE, nrow = 20, ncol = 5))
    expect_identical(as.matrix(s1), as.matrix(s2))

    s3 <- simil(mat1_test, mat2_test, margin = 2,
                mask = Matrix(1.0, nrow = 20, ncol = 5))
    expect_identical(as.matrix(s1), as.matrix(s3))
})
