// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// cpp_linear
S4 cpp_linear(arma::sp_mat& mt1, arma::sp_mat& mt2, arma::sp_mat& mask, const int method, unsigned int rank, const double limit, bool symm, const bool drop0, const bool use_nan, const bool use_mask, const bool sparse, const int digits, const int thread);
RcppExport SEXP _proxyC_cpp_linear(SEXP mt1SEXP, SEXP mt2SEXP, SEXP maskSEXP, SEXP methodSEXP, SEXP rankSEXP, SEXP limitSEXP, SEXP symmSEXP, SEXP drop0SEXP, SEXP use_nanSEXP, SEXP use_maskSEXP, SEXP sparseSEXP, SEXP digitsSEXP, SEXP threadSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat& >::type mt1(mt1SEXP);
    Rcpp::traits::input_parameter< arma::sp_mat& >::type mt2(mt2SEXP);
    Rcpp::traits::input_parameter< arma::sp_mat& >::type mask(maskSEXP);
    Rcpp::traits::input_parameter< const int >::type method(methodSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type rank(rankSEXP);
    Rcpp::traits::input_parameter< const double >::type limit(limitSEXP);
    Rcpp::traits::input_parameter< bool >::type symm(symmSEXP);
    Rcpp::traits::input_parameter< const bool >::type drop0(drop0SEXP);
    Rcpp::traits::input_parameter< const bool >::type use_nan(use_nanSEXP);
    Rcpp::traits::input_parameter< const bool >::type use_mask(use_maskSEXP);
    Rcpp::traits::input_parameter< const bool >::type sparse(sparseSEXP);
    Rcpp::traits::input_parameter< const int >::type digits(digitsSEXP);
    Rcpp::traits::input_parameter< const int >::type thread(threadSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_linear(mt1, mt2, mask, method, rank, limit, symm, drop0, use_nan, use_mask, sparse, digits, thread));
    return rcpp_result_gen;
END_RCPP
}
// cpp_sd
NumericVector cpp_sd(arma::sp_mat& mt);
RcppExport SEXP _proxyC_cpp_sd(SEXP mtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat& >::type mt(mtSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_sd(mt));
    return rcpp_result_gen;
END_RCPP
}
// cpp_nz
NumericVector cpp_nz(arma::sp_mat& mt);
RcppExport SEXP _proxyC_cpp_nz(SEXP mtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat& >::type mt(mtSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_nz(mt));
    return rcpp_result_gen;
END_RCPP
}
// cpp_mask
S4 cpp_mask(IntegerVector v1, IntegerVector v2);
RcppExport SEXP _proxyC_cpp_mask(SEXP v1SEXP, SEXP v2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type v1(v1SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type v2(v2SEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_mask(v1, v2));
    return rcpp_result_gen;
END_RCPP
}
// cpp_pair
S4 cpp_pair(arma::sp_mat& mt1, arma::sp_mat& mt2, arma::sp_mat& mask, const int method, unsigned int rank, const double limit, const double weight, const double smooth, bool symm, const bool diag, const bool drop0, const bool use_nan, const bool use_mask, const bool sparse, const int digits, const int thread);
RcppExport SEXP _proxyC_cpp_pair(SEXP mt1SEXP, SEXP mt2SEXP, SEXP maskSEXP, SEXP methodSEXP, SEXP rankSEXP, SEXP limitSEXP, SEXP weightSEXP, SEXP smoothSEXP, SEXP symmSEXP, SEXP diagSEXP, SEXP drop0SEXP, SEXP use_nanSEXP, SEXP use_maskSEXP, SEXP sparseSEXP, SEXP digitsSEXP, SEXP threadSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat& >::type mt1(mt1SEXP);
    Rcpp::traits::input_parameter< arma::sp_mat& >::type mt2(mt2SEXP);
    Rcpp::traits::input_parameter< arma::sp_mat& >::type mask(maskSEXP);
    Rcpp::traits::input_parameter< const int >::type method(methodSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type rank(rankSEXP);
    Rcpp::traits::input_parameter< const double >::type limit(limitSEXP);
    Rcpp::traits::input_parameter< const double >::type weight(weightSEXP);
    Rcpp::traits::input_parameter< const double >::type smooth(smoothSEXP);
    Rcpp::traits::input_parameter< bool >::type symm(symmSEXP);
    Rcpp::traits::input_parameter< const bool >::type diag(diagSEXP);
    Rcpp::traits::input_parameter< const bool >::type drop0(drop0SEXP);
    Rcpp::traits::input_parameter< const bool >::type use_nan(use_nanSEXP);
    Rcpp::traits::input_parameter< const bool >::type use_mask(use_maskSEXP);
    Rcpp::traits::input_parameter< const bool >::type sparse(sparseSEXP);
    Rcpp::traits::input_parameter< const int >::type digits(digitsSEXP);
    Rcpp::traits::input_parameter< const int >::type thread(threadSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_pair(mt1, mt2, mask, method, rank, limit, weight, smooth, symm, diag, drop0, use_nan, use_mask, sparse, digits, thread));
    return rcpp_result_gen;
END_RCPP
}
// cpp_get_max_thread
int cpp_get_max_thread();
RcppExport SEXP _proxyC_cpp_get_max_thread() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(cpp_get_max_thread());
    return rcpp_result_gen;
END_RCPP
}
// cpp_tbb_enabled
bool cpp_tbb_enabled();
RcppExport SEXP _proxyC_cpp_tbb_enabled() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(cpp_tbb_enabled());
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_proxyC_cpp_linear", (DL_FUNC) &_proxyC_cpp_linear, 13},
    {"_proxyC_cpp_sd", (DL_FUNC) &_proxyC_cpp_sd, 1},
    {"_proxyC_cpp_nz", (DL_FUNC) &_proxyC_cpp_nz, 1},
    {"_proxyC_cpp_mask", (DL_FUNC) &_proxyC_cpp_mask, 2},
    {"_proxyC_cpp_pair", (DL_FUNC) &_proxyC_cpp_pair, 16},
    {"_proxyC_cpp_get_max_thread", (DL_FUNC) &_proxyC_cpp_get_max_thread, 0},
    {"_proxyC_cpp_tbb_enabled", (DL_FUNC) &_proxyC_cpp_tbb_enabled, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_proxyC(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
