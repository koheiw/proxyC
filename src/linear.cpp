#include "proxyc.h"
#include "dev.h"
using namespace proxyc;
using namespace arma;

rowvec nnz(const sp_mat& mt) {
    rowvec v(mt.n_cols, fill::zeros);
    if (mt.is_empty()) return(v);
    for (uword i = 0; i < mt.n_cols; i++) {
        v[i] = sum(colvec(mt.col(i)) != 0);
    }
    return(v);
}

rowvec stddev(const sp_mat& mt, const int norm_type) {
    rowvec v(mt.n_cols, fill::zeros);
    if (mt.is_empty()) return(v);
    for (uword i = 0; i < mt.n_cols; i++) {
        v[i] = stddev(colvec(mt.col(i)), norm_type);
    }
    return(v);
}

rowvec mean(const sp_mat& mt) {
    rowvec v(mt.n_cols, fill::zeros);
    if (mt.is_empty()) return(v);
    for (uword i = 0; i < mt.n_cols; i++) {
        v[i] = mean(colvec(mt.col(i)));
    }
    return(v);
}

void proxy_linear(const uword i,
                  const sp_mat& mt1t, const sp_mat& mt2, Triplets& simil_tri,
                  const rowvec& square1, const rowvec& center1,
                  const rowvec& square2, const rowvec& center2,
                  const int method,
                  const unsigned int rank, const double limit,
                  const bool symm, const bool drop0, const bool use_nan) {

    uword nrow = mt1t.n_rows;
    uword ncol = mt1t.n_cols;

    rowvec v1, v2;
    std::vector<double> simils(nrow);
    //for (uword i = begin; i < end; i++) {
        switch (method) {
        case 1: // cosine similarity
            simils = to_vector(trans(mt1t * mt2.col(i)) / (square1 * square2[i]));
            break;
        case 2: // correlation similarity
            v1 = rowvec(trans(mt1t * mt2.col(i)));
            v2 = center1 * center2[i] * ncol;
            simils = to_vector(((v1 - v2) / ncol) / (square1 * square2[i]));
            simils = replace_inf(simils);
            break;
        case 3: // euclidean distance
            simils = to_vector(sqrt(trans(mt1t * mt2.col(i)) * -2 + square1 + square2[i]));
            break;
        }
        double l = get_limit(simils, rank, limit);
        for (std::size_t k = 0; k < simils.size(); k++) {
            if (symm && k > i) continue;
            if (drop0 && simils[k] == 0) continue;
            if (simils[k] >= l || (use_nan && std::isnan(simils[k])))
                simil_tri.push_back(std::make_tuple(k, i, simils[k]));
        }
    //}
}

// [[Rcpp::export]]
S4 cpp_linear(arma::sp_mat& mt1,
              arma::sp_mat& mt2,
              const int method,
              unsigned int rank,
              const double limit = -1.0,
              bool symm = false,
              const bool drop0 = false,
              const bool use_nan = false,
              const int thread = -1) {

    if (mt1.n_rows != mt2.n_rows)
        throw std::range_error("Invalid matrix objects");

    uword ncol1 = mt1.n_cols;
    uword ncol2 = mt2.n_cols;
    if (rank < 1) rank = 1;
    symm = symm && rank == ncol2;

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
    std::size_t I = ncol2;
#if PROXYC_USE_TBB
    tbb::task_arena arena(thread);
    arena.execute([&]{
        tbb::parallel_for(tbb::blocked_range<int>(0, I), [&](tbb::blocked_range<int> r) {
            for (int i = r.begin(); i < r.end(); i++) {
                proxy_linear(i, mt1, mt2, simil_tri,
                             square1, center1, square2, center2,
                             method, rank, limit, symm, drop0, use_nan);
            }
        });
    });
# else
    for (std::size_t i = 0; i < I; i++) {
        proxy_linear(i, mt1, mt2, simil_tri,
                     square1, center1, square2, center2,
                     method, rank, limit, symm, drop0, use_nan);
    }
# endif

    return to_matrix(simil_tri, ncol1, ncol2, symm);
}

// [[Rcpp::export]]
NumericVector cpp_sd(arma::sp_mat& mt) {
    std::vector<double> sds = to_vector(stddev(mt, 0));
    return wrap(sds);
}

// [[Rcpp::export]]
NumericVector cpp_nz(arma::sp_mat& mt) {
    std::vector<double> nzs = to_vector(nnz(mt));
    return wrap(nzs);
}

/***R
mt <- Matrix::rsparsematrix(100, 100, 0.01)
system.time(
    out <- cpp_linear(mt, mt, 1, 100, symm = TRUE)
)
*/

