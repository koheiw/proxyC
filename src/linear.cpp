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
                  const sp_mat& mt1t, const sp_mat& mt2, const sp_mat& mask,
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

// [[Rcpp::export]]
S4 cpp_linear(arma::sp_mat& mt1,
              arma::sp_mat& mt2,
              const int method,
              arma::sp_mat& mask,
              unsigned int rank,
              const double limit = -1.0,
              bool symm = false,
              const bool drop0 = false,
              const bool use_nan = false,
              const bool use_mask = false,
              const bool sparse = true,
              const int digits = 14,
              const int thread = -1) {

    if (mt1.n_rows != mt2.n_rows)
        throw std::range_error("Invalid matrix objects");

    if (use_mask & (mask.n_rows != mt1.n_cols || mask.n_cols != mt2.n_cols))
        throw std::range_error("Invalid mask object");

    uword ncol1 = mt1.n_cols;
    uword ncol2 = mt2.n_cols;
    if (rank < 1) rank = 1;
    symm = symm && rank == ncol2;

    //dev::Timer timer;
    //dev::start_timer("Compute magnitude", timer);
    rowvec square1(ncol1), center1(ncol1), sum1(ncol1);
    rowvec square2(ncol2), center2(ncol2), sum2(ncol2);
    switch (method) {
    case 2: // cosine
        square1 = rowvec(sqrt(mat(sum(mt1 % mt1, 0))));
        square2 = rowvec(sqrt(mat(sum(mt2 % mt2, 0))));
        break;
    case 3: // correlation
        square1 = stddev(mt1, 1);
        center1 = mean(mt1);
        square2 = stddev(mt2, 1);
        center2 = mean(mt2);
        break;
    case 4: // euclidean distance
        square1 = rowvec(mat(sum(mt1 % mt1, 0)));
        square2 = rowvec(mat(sum(mt2 % mt2, 0)));
        break;
    case 5: // dice coefficient
        sum1 = sum(mt1, 0);
        sum2 = sum(mt2, 0);
        break;
    case 6: // edice coefficient
        sum1 = sum(square(mt1), 0);
        sum2 = sum(square(mt2), 0);
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
                proxy_linear(i, mt1, mt2, mask, simil_tri,
                             square1, center1, sum1, square2, center2, sum2,
                             method, rank, limit, symm, drop0, use_nan, use_mask, digits);
            }
        });
    });
# else
    for (std::size_t i = 0; i < I; i++) {
        proxy_linear(i, mt1, mt2, mask, simil_tri,
                     square1, center1, sum1, square2, center2, sum2,
                     method, rank, limit, symm, drop0, use_nan, use_mask, digits);
    }
# endif

    return to_matrix(simil_tri, ncol1, ncol2, symm, sparse);
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
mt1 <- Matrix::rsparsematrix(100, 2, 1)
mt2 <- Matrix::rsparsematrix(100, 20, 1)
mask <- Matrix::rsparsematrix(2, 20, 0.1)
system.time(
    out <- cpp_linear(as(mt1 > 0, "dMatrix"), as(mt2 > 0, "dMatrix"), 4, mask, 100, symm = TRUE, use_mask = TRUE, drop0 = TRUE)
)
*/

