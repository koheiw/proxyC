require(Matrix)

set.seed(42)
smat_test1 <- rsparsematrix(100, 200, 0.9)
dmat_test1 <- as.matrix(smat_test1)

smat_test2 <- rsparsematrix(100, 200, 0.1)
dmat_test2 <- as.matrix(smat_test2)

test_that("simil the same results", {

    for (m in eval(formals(simil)$method)) {
        expect_equal(
            as.matrix(simil(smat_test1, method = m)),
            as.matrix(simil(dmat_test1, method = m))
        )

        expect_equal(
            as.matrix(simil(smat_test1, dmat_test1, method = m)),
            as.matrix(simil(dmat_test1, smat_test1, method = m))
        )

        expect_equal(
            as.matrix(simil(smat_test2, method = m)),
            as.matrix(simil(dmat_test2, method = m))
        )

        expect_equal(
            as.matrix(simil(smat_test2, dmat_test2, method = m)),
            as.matrix(simil(dmat_test2, smat_test2, method = m))
        )
    }

})

test_that("dist the same results", {

    for (m in eval(formals(dist)$method)) {
        expect_equal(
            as.matrix(dist(smat_test1, method = m)),
            as.matrix(dist(dmat_test1, method = m))
        )

        expect_equal(
            as.matrix(dist(smat_test1, dmat_test1, method = m)),
            as.matrix(dist(dmat_test1, smat_test1, method = m))
        )

        expect_equal(
            as.matrix(dist(smat_test2, method = m)),
            as.matrix(dist(dmat_test2, method = m))
        )

        expect_equal(
            as.matrix(dist(smat_test2, dmat_test2, method = m)),
            as.matrix(dist(dmat_test2, smat_test2, method = m))
        )
    }

})

# positive values keep fjaccard in [0, 1]; items are columns, transposed for margin = 1
sx <- rsparsematrix(60, 40, 0.5, rand.x = function(n) runif(n, 0.1, 1))
sy <- rsparsematrix(60, 30, 0.5, rand.x = function(n) runif(n, 0.1, 1))
dimnames(sx) <- list(paste0("f", 1:60), paste0("x", 1:40))
dimnames(sy) <- list(paste0("f", 1:60), paste0("y", 1:30))
inputs <- function(margin) {
    if (margin == 1)
        list(sx = t(sx), sy = t(sy), dx = as.matrix(t(sx)), dy = as.matrix(t(sy)))
    else
        list(sx = sx, sy = sy, dx = as.matrix(sx), dy = as.matrix(sy))
}

test_that("dense and sparse input give the same result for both margins", {
    for (fun in list(simil, dist)) {
        for (m in eval(formals(fun)$method)) {
            for (margin in 1:2) {
                i <- inputs(margin)
                expected <- as.matrix(fun(i$sx, i$sy, margin = margin, method = m))
                expect_equal(as.matrix(fun(i$dx, i$dy, margin = margin, method = m)),
                             expected, tolerance = 1e-6, info = paste(m, margin))
                expect_equal(as.matrix(fun(i$dx, i$sy, margin = margin, method = m)),
                             expected, tolerance = 1e-6, info = paste(m, margin, "mixed"))
                expect_equal(as.matrix(fun(i$dx, margin = margin, method = m)),
                             as.matrix(fun(i$sx, margin = margin, method = m)),
                             tolerance = 1e-6, info = paste(m, margin, "y = NULL"))
            }
        }
    }
})

test_that("dimnames follow the margin for dense input", {
    i <- inputs(1)
    s <- simil(i$dx, i$dy, margin = 1)
    expect_identical(dimnames(s), list(rownames(i$dx), rownames(i$dy)))
    i <- inputs(2)
    s <- simil(i$dx, i$dy, margin = 2)
    expect_identical(dimnames(s), list(colnames(i$dx), colnames(i$dy)))
})

test_that("diag, rank, min_simil and mask work with dense input", {
    for (margin in 1:2) {
        i <- inputs(margin)
        expect_equal(as.matrix(simil(i$dx, i$dx, margin = margin, method = "jaccard", diag = TRUE)),
                     as.matrix(simil(i$sx, i$sx, margin = margin, method = "jaccard", diag = TRUE)))
        expect_equal(as.matrix(simil(i$dx, i$dy, margin = margin, rank = 3)),
                     as.matrix(simil(i$sx, i$sy, margin = margin, rank = 3)))
        expect_equal(as.matrix(simil(i$dx, i$dy, margin = margin, min_simil = 0.6)),
                     as.matrix(simil(i$sx, i$sy, margin = margin, min_simil = 0.6)))
        msk <- mask(rep(letters[1:4], 10), rep(letters[1:3], 10))
        expect_equal(as.matrix(simil(i$dx, i$dy, margin = margin, mask = msk)),
                     as.matrix(simil(i$sx, i$sy, margin = margin, mask = msk)))
        expect_error(simil(i$dx, i$dy, margin = margin, mask = t(msk)),
                     "The shape of mask must be 40 x 30")
    }
})

test_that("integer, logical and Matrix dense input give the same result as double input", {
    d <- as.matrix(inputs(2)$sx > 0.5)
    expected <- as.matrix(simil(d * 1, margin = 2, method = "jaccard"))
    expect_equal(as.matrix(simil(d, margin = 2, method = "jaccard")), expected)
    storage.mode(d) <- "integer"
    expect_equal(as.matrix(simil(d, margin = 2, method = "jaccard")), expected)
    expect_equal(as.matrix(simil(Matrix(d * 1, sparse = FALSE), margin = 2, method = "jaccard")), expected)
})

test_that("dense input is not modified", {
    i <- inputs(1)
    before <- i$dx + 0
    simil(i$dx, margin = 1, method = "jaccard")
    simil(i$dx, margin = 1, method = "cosine")
    expect_identical(i$dx, before)

    # a square matrix with margin = 2 is read in place, so a transpose there would hit it
    sq <- as.matrix(rsparsematrix(30, 30, 0.5))
    before <- sq + 0
    simil(sq, margin = 2, method = "cosine")
    dist(sq, margin = 2, method = "euclidean")
    expect_identical(sq, before)
})
