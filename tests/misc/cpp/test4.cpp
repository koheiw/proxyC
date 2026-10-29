#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
void test(arma::sp_mat& x) {
    arma::sp_mat y = arma::trans(x);
}

// [[Rcpp::export]]
void test_dense(arma::mat& x) {
    arma::mat y = arma::trans(x);
}

// [[Rcpp::export]]
void test_dense2(NumericMatrix& x) {
    NumericMatrix y = transpose(x);
}

/*** R
smat <- Matrix::rsparsematrix(1000, 100000, density = 0.5)
dmat <- as.matrix(smat)

microbenchmark::microbenchmark(
    test(smat),
    test_dense(dmat),
    test_dense2(dmat),
    Matrix::t(smat),
    t(dmat),
    times = 10
)
*/
