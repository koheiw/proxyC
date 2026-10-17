#include "linear.h"
#include "proxyc.h"
#include "dev.h"
using namespace proxyc;
using namespace arma;

// [[Rcpp::export]]
S4 cpp_linear(SEXP s1,
              SEXP s2,
              const int method,
              arma::sp_mat& mask,
              unsigned int rank,
              const double limit = -1.0,
              bool symm = false,
              const bool drop0 = false,
              const bool use_nan = false,
              const bool use_mask = false,
              const bool sparse = true,
              const int digits = 14,
              const int thread = -1) {

    if (Rf_isS4(s1) && Rf_isS4(s2)) {
        arma::sp_mat mt1 = Rcpp::as<arma::sp_mat>(s1);
        arma::sp_mat mt2 = Rcpp::as<arma::sp_mat>(s2);
        return iter_linear(mt1, mt2, method, mask, rank, limit, symm, drop0,
                           use_nan, use_mask, sparse, digits, thread);
    } else if (Rf_isS4(s1) && !Rf_isS4(s2)) {
        arma::sp_mat mt1 = Rcpp::as<arma::sp_mat>(s1);
        arma::mat mt2 = Rcpp::as<arma::mat>(s2);
        return iter_linear(mt1, mt2, method, mask, rank, limit, symm, drop0,
                           use_nan, use_mask, sparse, digits, thread);
    } else if (!Rf_isS4(s1) && Rf_isS4(s2)) {
        arma::mat mt1 = Rcpp::as<arma::mat>(s1);
        arma::sp_mat mt2 = Rcpp::as<arma::sp_mat>(s2);
        return iter_linear(mt1, mt2, method, mask, rank, limit, symm, drop0,
                           use_nan, use_mask, sparse, digits, thread);
    } else {
        arma::mat mt1 = Rcpp::as<arma::mat>(s1);
        arma::mat mt2 = Rcpp::as<arma::mat>(s2);
        return iter_linear(mt1, mt2, method, mask, rank, limit, symm, drop0,
                           use_nan, use_mask, sparse, digits, thread);
    }
}

// [[Rcpp::export]]
NumericVector cpp_sd(SEXP s) {
    if (Rf_isS4(s)) {
        arma::sp_mat mt = Rcpp::as<arma::sp_mat>(s);
        std::vector<double> sds = to_vector(stddev(mt, 0));
        return wrap(sds);
    } else {
        arma::mat mt = Rcpp::as<arma::mat>(s);
        std::vector<double> sds = to_vector(stddev(mt, 0));
        return wrap(sds);
    }
}

// [[Rcpp::export]]
NumericVector cpp_nz(SEXP s) {
    if (Rf_isS4(s)) {
        arma::sp_mat mt = Rcpp::as<arma::sp_mat>(s);
        std::vector<double> nzs = to_vector(nnz(mt));
        return wrap(nzs);
    } else {
        arma::mat mt = Rcpp::as<arma::mat>(s);
        std::vector<double> nzs = to_vector(nnz(mt));
        return wrap(nzs);
    }
}


/***R
mt1 <- Matrix::rsparsematrix(100, 2, 1)
mt2 <- Matrix::rsparsematrix(100, 20, 1)
mask <- Matrix::rsparsematrix(2, 20, 0.1)
cpp_linear(mt1, mt2, 4, mask, 100, symm = TRUE, use_mask = TRUE, drop0 = TRUE)
cpp_linear(as.matrix(mt1), as.matrix(mt2), 4, mask, 100, symm = TRUE, use_mask = TRUE, drop0 = TRUE)

*/

