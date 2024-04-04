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

    inline S4 to_matrix(Triplets& tri, int nrow, int ncol, bool symmetric) {

        std::size_t l = tri.size();
        IntegerVector dim_ = IntegerVector::create(nrow, ncol);
        IntegerVector i_(l), j_(l);
        NumericVector x_(l);

        for (std::size_t k = 0; k < tri.size(); k++) {
            i_[k] = std::get<0>(tri[k]);
            j_[k] = std::get<1>(tri[k]);
            x_[k] = std::get<2>(tri[k]);
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

    inline double get_limit(std::vector<double> simils, const unsigned int rank, double limit) {

        if (simils.size() > rank) {
            std::nth_element(simils.begin(), simils.begin() + rank - 1, simils.end(),
                             std::greater<double>());
            if (limit < simils[rank - 1])
                limit = simils[rank - 1];
        }
        return limit;
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
