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
              const int thread = -1,
              const int margin = 2,
              const bool boolean = false) {

    return with_items(s1, margin, boolean, [&](const auto& mt1) {
        return with_items(s2, margin, boolean, [&](const auto& mt2) {
            return iter_linear(mt1, mt2, method, mask, rank, limit, symm, drop0,
                               use_nan, use_mask, sparse, digits, thread);
        });
    });
}

// [[Rcpp::export]]
NumericVector cpp_sd(SEXP s, const int margin = 2) {
    return with_items(s, margin, false, [](const auto& mt) {
        return NumericVector(wrap(to_vector(stddev(mt, 0))));
    });
}

// [[Rcpp::export]]
NumericVector cpp_nz(SEXP s, const int margin = 2) {
    return with_items(s, margin, false, [](const auto& mt) {
        return NumericVector(wrap(to_vector(nnz(mt))));
    });
}


/***R
mt1 <- Matrix::rsparsematrix(100, 2, 1)
mt2 <- Matrix::rsparsematrix(100, 20, 1)
mask <- Matrix::rsparsematrix(2, 20, 0.1)
cpp_linear(mt1, mt2, 4, mask, 100, symm = TRUE, use_mask = TRUE, drop0 = TRUE)
cpp_linear(as.matrix(mt1), as.matrix(mt2), 4, mask, 100, symm = TRUE, use_mask = TRUE, drop0 = TRUE)

*/

