#include "proxyc.h"
#include "utility.h"

// [[Rcpp::export]]
int cpp_get_max_thread() {
#if PROXYC_USE_TBB
    return tbb::this_task_arena::max_concurrency();
#else
    return 1;
#endif
}

// [[Rcpp::export]]
bool cpp_tbb_enabled(){
#if PROXYC_USE_TBB
    return true;
#else
    return false;
#endif
}

// [[Rcpp::export]]
IntegerVector cpp_nz(SEXP s, const int margin = 2) {
    if (Rf_isS4(s)) {
        const arma::sp_mat mt = to_sparse(s, margin, false);
        return wrap(to_vector(nnz(mt)));
    } else {
        const arma::mat mt = to_dense(s, margin, false);
        return wrap(to_vector(nnz(mt)));
    }
}

// [[Rcpp::export]]
NumericVector cpp_sd(SEXP s, const bool use_nan = false, const int margin = 2) {
    if (Rf_isS4(s)) {
        const arma::sp_mat mt = to_sparse(s, margin, false);
        return wrap(to_vector(stddev(mt, 0, use_nan)));
    } else {
        const arma::mat mt = to_dense(s, margin, false);
        return wrap(to_vector(stddev(mt, 0)));
    }
}

// [[Rcpp::export]]
NumericVector cpp_max(SEXP s, const bool use_nan = false, const int margin = 2) {
    if (Rf_isS4(s)) {
        const arma::sp_mat mt = to_sparse(s, margin, false);
        return wrap(to_vector(max(mt, use_nan)));
    } else {
        const arma::mat mt = to_dense(s, margin, false);
        return wrap(to_vector(max(mt)));
    }
}

// [[Rcpp::export]]
NumericVector cpp_min(SEXP s, const bool use_nan = false, const int margin = 2) {
    if (Rf_isS4(s)) {
        const arma::sp_mat mt = to_sparse(s, margin, false);
        return wrap(to_vector(min(mt, use_nan)));
    } else {
        const arma::mat mt = to_dense(s, margin, false);
        return wrap(to_vector(min(mt)));
    }
}

// [[Rcpp::export]]
IntegerVector cpp_maxrow(SEXP s, const int mode = 1,
                         const bool use_nan = false, const int margin = 2) {
    if (Rf_isS4(s)) {
        const arma::sp_mat mt = to_sparse(s, margin, false);
        return wrap(to_vector(maxrow(mt, mode, use_nan)));
    } else {
        const arma::mat mt = to_dense(s, margin, false);
        return wrap(to_vector(maxrow(mt, mode)));
    }
}
