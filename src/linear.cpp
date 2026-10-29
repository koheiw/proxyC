#include "linear.h"
#include "proxyc.h"
#include "dev.h"
using namespace proxyc;
using namespace arma;

// [[Rcpp::export]]
S4 cpp_linear(const SEXP s1,
              const SEXP s2,
              const int method,
              const arma::sp_mat& mask,
              const unsigned int rank,
              const double limit = -1.0,
              const bool symm = false,
              const bool drop0 = false,
              const bool use_nan = false,
              const bool use_mask = false,
              const bool sparse = true,
              const int digits = 14,
              const int thread = -1) {

    if (Rf_isS4(s1) && Rf_isS4(s2)) {
        arma::sp_mat mt1 = to_sparse(s1);
        arma::sp_mat mt2 = to_sparse(s2);
        return iter_linear(mt1, mt2, method, mask, rank, limit, symm, drop0,
                           use_nan, use_mask, sparse, digits, thread);
    } else if (Rf_isS4(s1) && !Rf_isS4(s2)) {
        arma::sp_mat mt1 = to_sparse(s1);
        arma::mat mt2 = to_dense(s2);
        return iter_linear(mt1, mt2, method, mask, rank, limit, symm, drop0,
                           use_nan, use_mask, sparse, digits, thread);
    } else if (!Rf_isS4(s1) && Rf_isS4(s2)) {
        arma::mat mt1 = to_dense(s1);
        arma::sp_mat mt2 = to_sparse(s2);
        return iter_linear(mt1, mt2, method, mask, rank, limit, symm, drop0,
                           use_nan, use_mask, sparse, digits, thread);
    } else {
        arma::mat mt1 = to_dense(s1);
        arma::mat mt2 = to_dense(s2);
        return iter_linear(mt1, mt2, method, mask, rank, limit, symm, drop0,
                           use_nan, use_mask, sparse, digits, thread);
    }
}


/***R
mt1 <- Matrix::rsparsematrix(100, 2, 1)
mt2 <- Matrix::rsparsematrix(100, 20, 1)
mask <- Matrix::rsparsematrix(2, 20, 0.1)
cpp_linear(mt1, mt2, 4, mask, 100, symm = TRUE, use_mask = TRUE, drop0 = TRUE)
cpp_linear(as.matrix(mt1), as.matrix(mt2), 4, mask, 100, symm = TRUE, use_mask = TRUE, drop0 = TRUE)

*/

