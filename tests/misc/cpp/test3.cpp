// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
void test_colvec(arma::sp_mat& x) {
    for (uword i = 0; i < x.n_cols; i++) {
        const colvec v = colvec(x.col(i));
    }
}

// [[Rcpp::export]]
void test_as_dense(arma::sp_mat& x) {
    for (uword i = 0; i < x.n_cols; i++) {
        const colvec v = x.col(i).as_dense();
    }
}


/*** R
mat <- Matrix::rsparsematrix(100000, 1000, density = 0.1)
microbenchmark::microbenchmark(
    test_colvec(mat),
    test_as_dense(mat),
    times = 100
)
*/
