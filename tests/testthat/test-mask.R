require(Matrix)

set.seed(1234)
smat1_test <- rsparsematrix(50, 20, 0.5)
colnames(smat1_test) <- rep(c("a", "b", "c", "d", "e"), each = 4)
smat2_test <- rsparsematrix(50, 5, 0.5)
colnames(smat2_test) <- c("a", "b", "c", "d", "e")
mmat_test <- mask(colnames(smat1_test), colnames(smat2_test))

test_that("mask is working", {

    c1 <- rep(c("a", "b", "c", "d", "e"), each = 4)
    c2 <- c("a", "b", "c", "d", "e")

    mmat1 <- mask(c1, c2)
    expect_equivalent(class(mmat1), "lgTMatrix")
    expect_identical(rownames(mmat1), c1)
    expect_identical(colnames(mmat1), c2)

    expect_error(
        mask(c1, as.factor(c2)),
        "x and y must be the same type of vectors"
    )

    n1 <- rep(1:5, each = 4)
    n2 <- 1:5

    mmat2 <- mask(n1, n2)
    expect_equivalent(class(mmat2), "lgTMatrix")
    expect_null(rownames(mmat2))
    expect_null(colnames(mmat2))

    expect_error(
        mask(n1, as.character(n2)),
        "x and y must be the same type of vectors"
    )
})

test_that("mask works with simil linear", {

    # masked dense
    s1 <- simil(smat1_test, smat2_test, margin = 2,
                mask = mmat_test)
    expect_identical(as.matrix(mmat_test != 0), as.matrix(s1 != 0))

    # masked ranked
    s2 <- simil(smat1_test, smat2_test, margin = 2, rank = 1, drop0 = TRUE,
                mask = mmat_test)
    expect_identical(colSums(s2), sapply(split(s2@x, rownames(s2)[s2@i + 1]), max))

    # masked min
    s3 <- simil(smat1_test, smat2_test, margin = 2, min_simil = 0, drop0 = TRUE,
                mask = mmat_test)
    expect_true(all(sapply(split(s3@x, rownames(s3)[s3@i + 1]), min) > 0))

    # masked min
    s4 <- simil(smat1_test, smat2_test, margin = 2, min_simil = -0.1, drop0 = TRUE,
                mask = mmat_test)
    expect_true(all(sapply(split(s4@x, rownames(s4)[s4@i + 1]), min) > -0.1))

    expect_error(
        simil(smat1_test, smat2_test, margin = 2, mask = mmat_test[,-1]),
        "The shape of mask must be 20 x 5"
    )
})

test_that("mask works with simil pair", {

    # masked dense
    s1 <- simil(smat1_test, smat2_test, margin = 2, method = "jaccard",
                mask = mmat_test)
    expect_identical(as.matrix(mmat_test != 0), as.matrix(s1 != 0))

    # masked ranked
    s2 <- simil(smat1_test, smat2_test, margin = 2, method = "jaccard", rank = 1, drop0 = TRUE,
                mask = mmat_test)
    expect_identical(colSums(s2), sapply(split(s2@x, rownames(s2)[s2@i + 1]), max))

    # masked min
    s3 <- simil(smat1_test, smat2_test, margin = 2,method = "jaccard",  min_simil = 0, drop0 = TRUE,
                mask = mmat_test)
    expect_true(all(sapply(split(s3@x, rownames(s3)[s3@i + 1]), min) > 0))

    # masked min
    s4 <- simil(smat1_test, smat2_test, margin = 2, method = "jaccard", min_simil = -0.1, drop0 = TRUE,
                mask = mmat_test)
    expect_true(all(sapply(split(s4@x, rownames(s4)[s4@i + 1]), min) > -0.1))

    expect_error(
        simil(smat1_test, smat2_test, margin = 2, method = "jaccard", mask = mmat_test[,-1]),
        "The shape of mask must be 20 x 5"
    )
})

