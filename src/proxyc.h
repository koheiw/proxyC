#ifndef PROXYC // prevent redefining
#define PROXYC

#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>
#ifdef TBB
#include <tbb/tbb.h>
#ifdef ONETBB_SPEC_VERSION
using namespace oneapi; // only Windows R 4.3.x or later
#endif
#define PROXYC_USE_TBB true
#else
#define PROXYC_USE_TBB false
#endif

// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;
using namespace std;

// setting for unordered_map and unordered_set
const float GLOBAL_PATTERNS_MAX_LOAD_FACTOR = 0.1;
const float GLOBAL_NGRAMS_MAX_LOAD_FACTOR = 0.5;

namespace proxyc{

#if PROXYC_USE_TBB
    typedef std::tuple<unsigned int, unsigned int, double> Triplet;
    typedef tbb::concurrent_vector<Triplet> Triplets;
#else
    typedef std::tuple<unsigned int, unsigned int, double> Triplet;
    typedef std::vector<Triplet> Triplets;
#endif

    inline S4 to_matrix(Triplets& tri, int nrow, int ncol, bool symmetric, bool sparse) {

        IntegerVector dim_ = IntegerVector::create(nrow, ncol);
        if (!sparse) {
            if (symmetric) {
                std::size_t l = nrow * (nrow + 1) / 2;
                NumericVector x_(l, 0);
                for (Triplet t : tri) {
                    std::size_t k = std::get<0>(t) + (std::get<1>(t) * (std::get<1>(t) + 1) / 2);
                    x_[k] = std::get<2>(t);
                }
                S4 simil_("dspMatrix");
                simil_.slot("x") = x_;
                simil_.slot("Dim") = dim_;
                simil_.slot("uplo") = "U";
                return simil_;
            } else {
                std::size_t l = nrow * ncol;
                NumericVector x_(l, 0);
                for (Triplet t : tri) {
                    std::size_t k = std::get<0>(t) + (std::get<1>(t) * nrow);
                    x_[k] = std::get<2>(t);
                }
                S4 simil_("dgeMatrix");
                simil_.slot("x") = x_;
                simil_.slot("Dim") = dim_;
                return simil_;
            }
        } else {
            std::size_t l = tri.size();
            NumericVector x_(l, 0);
            IntegerVector i_(l), j_(l);
            std::size_t k = 0;
            for (Triplet t : tri) {
                i_[k] = std::get<0>(t);
                j_[k] = std::get<1>(t);
                x_[k] = std::get<2>(t);
                k++;
            }
            if (symmetric) {
                S4 simil_("dsTMatrix");
                simil_.slot("i") = i_;
                simil_.slot("j") = j_;
                simil_.slot("x") = x_;
                simil_.slot("Dim") = dim_;
                simil_.slot("uplo") = "U";
                return simil_;

            } else {
                S4 simil_("dgTMatrix");
                simil_.slot("i") = i_;
                simil_.slot("j") = j_;
                simil_.slot("x") = x_;
                simil_.slot("Dim") = dim_;
                return simil_;
            }
        }
    }

    inline double get_limit(std::vector<double> simils, const unsigned int rank,
                            double limit) {

        if (simils.size() > rank) {
            std::sort(simils.begin(), simils.end(),
                      [] (auto x, auto y) {
                          // move nan to the last
                          if (std::isnan(x)) return false;
                          if (std::isnan(y)) return true;
                          return x > y;
                      });
            double s = simils[rank - 1];
            if (limit < s && !std::isnan(s))
                limit = s;
        }
        return limit;
    }

    inline std::vector<double> round(std::vector<double> simils, const int digits) {
        double shift = std::pow(10, digits);
        for (auto it = simils.begin() ; it != simils.end(); ++it) {
            *it = std::round(*it * shift) / shift;
        }
        return simils;
    }

    inline std::vector<double> replace_masked(std::vector<double> simils,
                                           arma::sp_vec mask) {
        std::size_t j = 0;
        for (auto it = simils.begin() ; it != simils.end(); ++it) {
            if (mask[j] == 0)
                *it = std::numeric_limits<double>::quiet_NaN();
            j++;
        }
        return simils;
    }

    inline std::vector<double> replace_nan(std::vector<double> simils) {
        for (auto it = simils.begin() ; it != simils.end(); ++it) {
            if (std::isnan(*it)) {
                *it = 0;
            }
        }
        return simils;
    }

    inline std::vector<double> replace_inf(std::vector<double> simils) {
        for (auto it = simils.begin() ; it != simils.end(); ++it) {
            if (std::isinf(*it)) {
                *it = std::numeric_limits<double>::quiet_NaN();
            }
        }
        return simils;
    }

    inline double replace_inf(double simil) {
        if (std::isinf(simil))
            return std::numeric_limits<double>::quiet_NaN();
        return simil;
    }

    inline std::vector<double> to_vector(const arma::sp_mat& mt) {
        return arma::conv_to< std::vector<double> >::from(arma::mat(mt));
    }

    inline std::vector<double> to_vector(const arma::rowvec& v) {
        return arma::conv_to< std::vector<double> >::from(v);
    }
}

#endif
