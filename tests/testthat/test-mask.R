require(Matrix)

set.seed(1234)
smat1 <- rsparsematrix(50, 20, 0.5)
colnames(smat1) <- rep(c("a", "b", "c", "d", "e"), each = 4)
smat2 <- rsparsematrix(50, 5, 0.5)
colnames(smat2) <- c("a", "b", "c", "d", "e")
mmat <- mask(colnames(smat1), colnames(smat2))

test_that("mask is working", {
    expect_equivalent(class(mmat), "lgTMatrix")
    expect_identical(rownames(mmat), colnames(smat1))
    expect_identical(colnames(mmat), colnames(smat2))
})

test_that("mask works with simil linear", {

    # masked dense
    s1 <- simil(smat1, smat2, margin = 2,
                        mask = mmat)
    expect_identical(as.matrix(mmat != 0), as.matrix(s1 != 0))

    # masked ranked
    s2 <- simil(smat1, smat2, margin = 2, rank = 1, drop0 = TRUE,
                        mask = mmat)
    expect_identical(colSums(s2), sapply(split(s2@x, rownames(s2)[s2@i + 1]), max))

    # masked min
    s3 <- simil(smat1, smat2, margin = 2, min_simil = 0, drop0 = TRUE,
                        mask = mmat)
    expect_true(all(sapply(split(s3@x, rownames(s3)[s3@i + 1]), min) > 0))

    # masked min
    s4 <- simil(smat1, smat2, margin = 2, min_simil = -0.1, drop0 = TRUE,
                        mask = mmat)
    expect_true(all(sapply(split(s4@x, rownames(s4)[s4@i + 1]), min) > -0.1))
})

test_that("mask works with simil pair", {

    # masked dense
    s1 <- simil(smat1, smat2, margin = 2, method = "jaccard",
                mask = mmat)
    expect_identical(as.matrix(mmat != 0), as.matrix(s1 != 0))

    # masked ranked
    s2 <- simil(smat1, smat2, margin = 2, method = "jaccard", rank = 1, drop0 = TRUE,
                        mask = mmat)
    expect_identical(colSums(s2), sapply(split(s2@x, rownames(s2)[s2@i + 1]), max))

    # masked min
    s3 <- simil(smat1, smat2, margin = 2,method = "jaccard",  min_simil = 0, drop0 = TRUE,
                        mask = mmat)
    expect_true(all(sapply(split(s3@x, rownames(s3)[s3@i + 1]), min) > 0))

    # masked min
    s4 <- simil(smat1, smat2, margin = 2, method = "jaccard", min_simil = -0.1, drop0 = TRUE,
                        mask = mmat)
    expect_true(all(sapply(split(s4@x, rownames(s4)[s4@i + 1]), min) > -0.1))
})

