require(Matrix)

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
