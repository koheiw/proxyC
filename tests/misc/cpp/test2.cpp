#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
void test_safe(arma::mat& x) {
    for (uword i = 0; i < x.n_cols; i++) {
        colvec v = x.col(i);
    }
}

// [[Rcpp::export]]
void test_safe_const(arma::mat& x) {
    for (uword i = 0; i < x.n_cols; i++) {
        const colvec v = x.col(i);
    }
}

// [[Rcpp::export]]
void test_unsafe(arma::mat& x) {
    for (uword i = 0; i < x.n_cols; i++) {
        colvec v = x.unsafe_col(i);
    }
}

// [[Rcpp::export]]
void test_unsafe_const(arma::mat& x) {
    for (uword i = 0; i < x.n_cols; i++) {
        const colvec v = x.unsafe_col(i);
    }
}


/*** R
mat <- as.matrix(Matrix::rsparsematrix(1000, 100000, density = 1))
microbenchmark::microbenchmark(
    test_safe(mat),
    test_safe_const(mat),
    test_unsafe(mat),
    test_unsafe_const(mat),
)
*/
