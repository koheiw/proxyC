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
NumericVector cpp_nz(SEXP s) {
    if (Rf_isS4(s)) {
        arma::sp_mat mt = Rcpp::as<arma::sp_mat>(s);
        std::vector<double> v = to_vector(nnz(mt));
        return wrap(v);
    } else {
        arma::mat mt = Rcpp::as<arma::mat>(s);
        std::vector<double> v = to_vector(nnz(mt));
        return wrap(v);
    }
}

// [[Rcpp::export]]
NumericVector cpp_sd(SEXP s, const bool use_nan = false) {
    if (Rf_isS4(s)) {
        arma::sp_mat mt = Rcpp::as<arma::sp_mat>(s);
        std::vector<double> v = to_vector(stddev(mt, 0, use_nan));
        return wrap(v);
    } else {
        arma::mat mt = Rcpp::as<arma::mat>(s);
        std::vector<double> v = to_vector(stddev(mt, 0));
        return wrap(v);
    }
}

// [[Rcpp::export]]
NumericVector cpp_max(SEXP s, const bool use_nan = false) {
    if (Rf_isS4(s)) {
        arma::sp_mat mt = Rcpp::as<arma::sp_mat>(s);
        std::vector<double> v = to_vector(max(mt, use_nan));
        return wrap(v);
    } else {
        arma::mat mt = Rcpp::as<arma::mat>(s);
        std::vector<double> v = to_vector(max(mt));
        return wrap(v);
    }
}

// [[Rcpp::export]]
NumericVector cpp_min(SEXP s, const bool use_nan = false) {
    if (Rf_isS4(s)) {
        arma::sp_mat mt = Rcpp::as<arma::sp_mat>(s);
        std::vector<double> v = to_vector(min(mt, use_nan));
        return wrap(v);
    } else {
        arma::mat mt = Rcpp::as<arma::mat>(s);
        std::vector<double> v = to_vector(min(mt));
        return wrap(v);
    }
}

// [[Rcpp::export]]
IntegerVector cpp_maxrow(SEXP s, const bool use_nan = false) {
    if (Rf_isS4(s)) {
        arma::sp_mat mt = Rcpp::as<arma::sp_mat>(s);
        std::vector<double> v = to_vector(maxrow(mt, use_nan));
        return wrap(v);
    } else {
        arma::mat mt = Rcpp::as<arma::mat>(s);
        std::vector<double> v = to_vector(maxrow(mt));
        return wrap(v);
    }
}
