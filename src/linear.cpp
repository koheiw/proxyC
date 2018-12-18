#include "armadillo.h"
#include "proxyc.h"
#include "dev.h"
using namespace proxyc;
using namespace arma;

rowvec stddev(const sp_mat& mt, const int norm_type) {
    rowvec v(mt.n_cols);
    for (uword i = 0; i < mt.n_cols; i++) {
        v[i] = stddev(colvec(mt.col(i)), norm_type);
    }
    return(v);
}

rowvec mean(const sp_mat& mt) {
    rowvec v(mt.n_cols);
    for (uword i = 0; i < mt.n_cols; i++) {
        v[i] = mean(colvec(mt.col(i)));
    }
    return(v);
}

struct proxy_linear : public Worker {

    const arma::sp_mat& mt1t; // input
    const arma::sp_mat& mt2; // input
    Triplets& simil_tri; // output
    const rowvec& square1;
    const rowvec& center1;
    const rowvec& square2;
    const rowvec& center2;
    const int method;
    const unsigned int rank;
    const double limit;
    const bool symm;

    proxy_linear(const sp_mat& mt1t_, const sp_mat& mt2_, Triplets& simil_tri_,
                      const rowvec& square1_, const rowvec& center1_,
                      const rowvec& square2_, const rowvec& center2_,
                      const int method_,
                      const unsigned int rank_, const double limit_, const bool symm_) :
        mt1t(mt1t_), mt2(mt2_), simil_tri(simil_tri_),
        square1(square1_), center1(center1_), square2(square2_), center2(center2_),
        method(method_), rank(rank_), limit(limit_), symm(symm_) {}

    void operator()(std::size_t begin, std::size_t end) {

        uword nrow = mt1t.n_rows;
        uword ncol = mt1t.n_cols;
        rowvec v1, v2;
        std::vector<double> simils(nrow);
        for (uword i = begin; i < end; i++) {
            switch (method) {
            case 1: // cosine similarity
                simils = to_vector(trans(mt1t * mt2.col(i)) / (square1 * square2[i]));
                break;
            case 2: // correlation similarity
                v1 = rowvec(trans(mt1t * mt2.col(i)));
                v2 = center1 * center2[i] * ncol;
                simils = to_vector(((v1 - v2) / ncol) / (square1 * square2[i]));
                break;
            case 3: // euclidean distance
                simils = to_vector(sqrt(trans(mt1t * mt2.col(i)) * -2 + square1 + square2[i]));
                break;
            }
            double l = get_limit(simils, rank, limit);
            for (std::size_t k = 0; k < simils.size(); k++) {
                if (symm && k > i) continue;
                if (simils[k] >= l) {
                    simil_tri.push_back(std::make_tuple(k, i, simils[k]));
                }
            }
        }
    }
};

// [[Rcpp::export]]
S4 cpp_linear(arma::sp_mat& mt1,
              arma::sp_mat& mt2,
              const int method,
              unsigned int rank,
              double limit = -1.0) {

    if (mt1.n_rows != mt2.n_rows)
        throw std::range_error("Invalid matrix objects");

    uword ncol1 = mt1.n_cols;
    uword ncol2 = mt2.n_cols;
    if (rank < 1) rank = 1;
    bool symm = rank == ncol1 && rank == ncol2;

    //dev::Timer timer;
    //dev::start_timer("Compute magnitude", timer);
    rowvec square1(ncol1), center1(ncol1), square2(ncol2), center2(ncol2);
    switch (method) {
    case 1: // cosine
        square1 = rowvec(sqrt(mat(sum(mt1 % mt1, 0))));
        square2 = rowvec(sqrt(mat(sum(mt2 % mt2, 0))));
        break;
    case 2: // correlation
        square1 = stddev(mt1, 1);
        center1 = mean(mt1);
        square2 = stddev(mt2, 1);
        center2 = mean(mt2);
        break;
    case 3: // euclidean distance
        square1 = rowvec(mat(sum(mt1 % mt1, 0)));
        square2 = rowvec(mat(sum(mt2 % mt2, 0)));
        break;
    }

    //dev::stop_timer("Compute magnitude", timer);
    //dev::start_timer("Compute similarity", timer);

    Triplets simil_tri;
    mt1 = trans(mt1);
    proxy_linear proxy_linear(mt1, mt2, simil_tri,
                              square1, center1, square2, center2,
                              method, rank, limit, symm);
    parallelFor(0, ncol2, proxy_linear);
    //dev::stop_timer("Compute similarity", timer);

    return to_matrix(simil_tri, ncol1, ncol2, symm);
}

/***R
mt <- Matrix::rsparsematrix(100, 100, 0.01)
system.time(
    out <- cpp_linear(mt, mt, 1, 100)
)
*/

