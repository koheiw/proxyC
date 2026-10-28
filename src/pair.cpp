#include "pair.h"
#include "proxyc.h"
#include "dev.h"
using namespace proxyc;
using namespace arma;

// [[Rcpp::export]]
S4 cpp_pair(SEXP s1,
            SEXP s2,
                  const int method,
                  const arma::sp_mat& mask,
                  unsigned int rank,
                  const double limit = -1.0,
                  const double weight = 1.0,
                  const double smooth = 0,
                  bool symm = false,
                  const bool diag = false,
                  const bool drop0 = false,
                  const bool use_nan = false,
                  const bool use_mask = false,
                  const bool sparse = true,
                  const int digits = 14,
                  const int thread = -1) {

    if (Rf_isS4(s1) && Rf_isS4(s2)) {
        arma::sp_mat mt1 = Rcpp::as<arma::sp_mat>(s1);
        arma::sp_mat mt2 = Rcpp::as<arma::sp_mat>(s2);
        return iter_pair(mt1, mt2, method, mask, rank, limit,
                         weight, smooth, symm, diag, drop0,
                         use_nan, use_mask, sparse, digits, thread);
    } else if (Rf_isS4(s1) && !Rf_isS4(s2)) {
        arma::sp_mat mt1 = Rcpp::as<arma::sp_mat>(s1);
        arma::mat mt2 = Rcpp::as<arma::mat>(s2);
        return iter_pair(mt1, mt2, method, mask, rank, limit,
                         weight, smooth, symm, diag, drop0,
                         use_nan, use_mask, sparse, digits, thread);
    } else if (!Rf_isS4(s1) && Rf_isS4(s2)) {
        arma::mat mt1 = Rcpp::as<arma::mat>(s1);
        arma::sp_mat mt2 = Rcpp::as<arma::sp_mat>(s2);
        return iter_pair(mt1, mt2, method, mask, rank, limit,
                         weight, smooth, symm, diag, drop0,
                         use_nan, use_mask, sparse, digits, thread);
    } else {
        arma::mat mt1 = Rcpp::as<arma::mat>(s1);
        arma::mat mt2 = Rcpp::as<arma::mat>(s2);
        return iter_pair(mt1, mt2, method, mask, rank, limit,
                         weight, smooth, symm, diag, drop0,
                         use_nan, use_mask, sparse, digits, thread);
    }

}

/***R
# mt <- Matrix::rsparsematrix(1000, 500, 0.01)
# microbenchmark::microbenchmark(
# out <- cpp_pair(mt, mt, 6, 100, symm = TRUE, smooth = 0.0)
# )

mt1 <- as.matrix(Matrix::rsparsematrix(100, 10, 1))
mt2 <- as.matrix(Matrix::rsparsematrix(100, 20, 1))
mask <- Matrix::rsparsematrix(10, 20, 0.1)
system.time(
    out <- cpp_pair(mt1, mt2, 1, mask, 100, symm = TRUE, use_mask = TRUE, drop0 = TRUE)
)

*/
