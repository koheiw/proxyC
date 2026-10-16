#include "proxyc.h"
#include "dev.h"
using namespace proxyc;
using namespace arma;

// template <typename T>
// rowvec nnz(const T& mt) {
//     rowvec v(mt.n_cols, fill::zeros);
//     if (mt.is_empty()) return(v);
//     for (uword i = 0; i < mt.n_cols; i++) {
//         v[i] = arma::sum(colvec(mt.col(i)) != 0);
//     }
//     return(v);
// }
//
// template <typename T>
// rowvec stddev(const T& mt, const int norm_type) {
//     rowvec v(mt.n_cols, fill::zeros);
//     if (mt.is_empty()) return(v);
//     for (uword i = 0; i < mt.n_cols; i++) {
//         v[i] = arma::stddev(colvec(mt.col(i)), norm_type);
//     }
//     return(v);
// }
//
// template <typename T>
// rowvec mean(const T& mt) {
//     rowvec v(mt.n_cols, fill::zeros);
//     if (mt.is_empty()) return(v);
//     for (uword i = 0; i < mt.n_cols; i++) {
//         v[i] = arma::mean(colvec(mt.col(i)));
//     }
//     return(v);
// }

template <typename T>
void proxy_linear(const uword i,
                  const T& mt1t, const T& mt2, const sp_mat& mask,
                  Triplets& simil_tri,
                  const rowvec& square1, const rowvec& center1, const rowvec& sum1,
                  const rowvec& square2, const rowvec& center2, const rowvec& sum2,
                  const int method,
                  const unsigned int rank, const double limit, const bool symm,
                  const bool drop0, const bool use_nan, const bool use_mask,
                  const int digits) {

    uword nrow = mt1t.n_rows;
    uword ncol = mt1t.n_cols;

    rowvec v1, v2;
    std::vector<double> simils(nrow);
    //for (uword i = begin; i < end; i++) {
    switch (method) {
    case 1: // cross-product
        simils = to_vector(trans(mt1t * mt2.col(i)));
        break;
    case 2: // cosine similarity
        simils = to_vector(trans(mt1t * mt2.col(i)) / (square1 * square2[i]));
        break;
    case 3: // correlation similarity
        v1 = rowvec(trans(mt1t * mt2.col(i)));
        v2 = center1 * center2[i] * ncol;
        simils = to_vector(((v1 - v2) / ncol) / (square1 * square2[i]));
        simils = replace_inf(simils);
        break;
    case 4: // euclidean distance
        simils = to_vector(sqrt(trans(mt1t * mt2.col(i)) * -2 + square1 + square2[i]));
        break;
    case 5: // dice coefficient
    case 6: // edice coefficient
        simils = to_vector(trans(mt1t * mt2.col(i) * 2) / (sum1 + sum2[i]));
        simils = replace_nan(simils);
        break;
    }
    if (use_mask)
        simils = replace_masked(simils, mask.col(i));
    simils = round(simils, digits);
    double l = get_limit(simils, rank, limit);
    for (std::size_t k = 0; k < simils.size(); k++) {
        double s = simils[k];
        if (symm && k > i) continue;
        if (drop0 && simils[k] == 0) continue;
        if (s >= l || (use_nan && std::isnan(s)))
            simil_tri.push_back(std::make_tuple(k, i, s));
    }
    //}
}
