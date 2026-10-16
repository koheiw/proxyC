#include "pair.h"

// [[Rcpp::export]]
S4 cpp_pair_sparse(arma::sp_mat& mt1,
                   arma::sp_mat& mt2,
                   const int method,
                   arma::sp_mat& mask,
                   unsigned int rank,
                   const double limit = -1.0,
                   const double weight = 1.0,
                   const double smooth = 0,
                   bool symm = false,
                   const bool diag = false,
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
    symm = symm && method != 12 && rank == ncol2; // exception for kullback

    Triplets simil_tri;
    std::size_t I = ncol2;
#if PROXYC_USE_TBB
    tbb::task_arena arena(thread);
    arena.execute([&]{
        tbb::parallel_for(tbb::blocked_range<int>(0, I), [&](tbb::blocked_range<int> r) {
            for (int i = r.begin(); i < r.end(); i++) {
                proxy_pair(i, mt1, mt2, mask, simil_tri, method, rank, limit, symm,
                           diag, weight, smooth, drop0, use_nan, use_mask, digits);
            }
        });
    });
#else
    for (std::size_t i = 0; i < I; i++) {
        proxy_pair(i, mt1, mt2, mask, simil_tri, method, rank, limit, symm,
                   diag, weight, smooth, drop0, use_nan, use_mask, digits);
    }
# endif

    return to_matrix(simil_tri, ncol1, ncol2, symm, sparse);

}

/***R
# mt <- Matrix::rsparsematrix(1000, 500, 0.01)
# microbenchmark::microbenchmark(
# out <- cpp_pair(mt, mt, 6, 100, symm = TRUE, smooth = 0.0)
# )

mt1 <- Matrix::rsparsematrix(100, 10, 1)
mt2 <- Matrix::rsparsematrix(100, 20, 1)
mask <- Matrix::rsparsematrix(10, 20, 0.1)
system.time(
    out <- cpp_pair_sparse(mt1, mt2, 1, mask, 100, symm = TRUE, use_mask = TRUE, drop0 = TRUE)
)

*/
