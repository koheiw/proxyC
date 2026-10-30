#include "proxyc.h"
#include "dev.h"
using namespace proxyc;
using namespace arma;

template <typename M1, typename M2>
void calc_linear(const uword i,
                 const M1& mt1t, const M2& mt2, const sp_mat& mask,
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
    double l = get_limit(simils, rank, limit); // NOTE: use arma::sort_index()?
    for (std::size_t k = 0; k < simils.size(); k++) {
        double s = simils[k];
        if (symm && k > i) continue;
        if (drop0 && simils[k] == 0) continue;
        if (s >= l || (use_nan && std::isnan(s)))
            simil_tri.push_back(std::make_tuple(k, i, s));
    }
    //}
}

template <typename M1, typename M2>
S4 iter_linear(M1& mt1,
                   const bool free_mt1,
                   const M2& mt2,
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
#if defined(__clang__) && __clang_major__ == 20
        // workaround for #65
        sum1 = sum(sp_mat(square(mt1)), 0);
        sum2 = sum(sp_mat(square(mt2)), 0);
#else
        sum1 = sum(square(mt1), 0);
        sum2 = sum(square(mt2), 0);
#endif
        break;
    }

    //dev::stop_timer("Compute magnitude", timer);
    //dev::start_timer("Compute similarity", timer);

    Triplets simil_tri;
    // only the transpose is used from here on
    M1 mt1t = trans(mt1);
    if (free_mt1)
        mt1.reset(); // not for a view on R's memory (see to_dense())
    std::size_t I = ncol2;
#if PROXYC_USE_TBB
    tbb::task_arena arena(thread);
    arena.execute([&]{
        tbb::parallel_for(tbb::blocked_range<int>(0, I), [&](tbb::blocked_range<int> r) {
            for (int i = r.begin(); i < r.end(); i++) {
                calc_linear(i, mt1t, mt2, mask, simil_tri,
                             square1, center1, sum1, square2, center2, sum2,
                             method, rank, limit, symm, drop0, use_nan, use_mask, digits);
            }
        });
    });
# else
    for (std::size_t i = 0; i < I; i++) {
        calc_linear(i, mt1t, mt2, mask, simil_tri,
                     square1, center1, sum1, square2, center2, sum2,
                     method, rank, limit, symm, drop0, use_nan, use_mask, digits);
    }
# endif

    return to_matrix(simil_tri, ncol1, ncol2, symm, sparse);
}

