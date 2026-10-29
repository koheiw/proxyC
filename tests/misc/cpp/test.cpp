// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
void test(arma::sp_mat& x) {
    Rcout << x.n_rows << " x "<< x.n_cols << "\n";
    x = arma::trans(x);
    Rcout << x.n_rows << " x "<< x.n_cols << "\n";
    x.reset();
    Rcout << x.n_rows << " x "<< x.n_cols << "\n";
}

// [[Rcpp::export]]
void test_dense(arma::mat& x) {
    Rcout << x.n_rows << " x "<< x.n_cols << "\n";
    x = arma::trans(x);
    Rcout << x.n_rows << " x "<< x.n_cols << "\n";
    x.reset();
    Rcout << x.n_rows << " x "<< x.n_cols << "\n";
}

// [[Rcpp::export]]
void test_dense2(NumericMatrix& x) {
    Rcout << x.nrow() << " x "<< x.ncol() << "\n";
    x = transpose(x);
    Rcout << x.nrow() << " x "<< x.ncol() << "\n";
    x = {};
    Rcout << x.nrow() << " x "<< x.ncol() << "\n";
}

/*** R
smat <- Matrix::rsparsematrix(10, 20, density = 0.5)
dmat <- as.matrix(smat)

dim(smat)
sum(smat != 0)
test(smat)
dim(smat)
sum(smat != 0)

dim(dmat)
sum(dmat != 0)
test_dense(dmat)
dim(dmat)
sum(dmat != 0)

dim(dmat)
sum(dmat != 0)
test_dense2(dmat)
dim(dmat)
sum(dmat != 0)
*/
