// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
void test1(const SEXP s) {
    arma::mat x = Rcpp::as<arma::mat>(s);
}

// [[Rcpp::export]]
void test2(const SEXP s) {
    NumericMatrix x_ = Rcpp::as<NumericMatrix>(s);
    arma::mat x(x_.begin(), x_.nrow(), x_.ncol(), false, true);
}

/*** R
smat <- Matrix::rsparsematrix(1000, 100000, density = 0.5)
dmat <- as.matrix(smat)

microbenchmark::microbenchmark(
    test1(dmat),
    test2(dmat),
    times = 10
)
*/
