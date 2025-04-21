#include "proxyc.h"
#include "dev.h"
using namespace proxyc;
using namespace arma;

// [[Rcpp::export]]
S4 cpp_mask(IntegerVector v1_,
            IntegerVector v2_,
            const int thread = -1) {

    std::vector<int> v1 = Rcpp::as<std::vector<int>>(v1_);
    std::vector<int> v2 = Rcpp::as<std::vector<int>>(v2_);

    std::size_t I = v1.size();
    std::size_t J = v2.size();
    Triplets simil_tri;

#if PROXYC_USE_TBB
    tbb::task_arena arena(thread);
    arena.execute([&]{
        tbb::parallel_for(tbb::blocked_range<int>(0, I), [&](tbb::blocked_range<int> r) {
            for (int i = r.begin(); i < r.end(); i++) {
                for (std::size_t j = 0; j < J; j++) {
                    if (v1[i] == v2[j])
                        simil_tri.push_back(std::make_tuple(i, j, 1.0));
                }
            }
        });
    });
# else
    for (std::size_t i = 0; i < I; i++) {
        for (std::size_t j = 0; j < J; j++) {
            if (v1[i] == v2[j])
                simil_tri.push_back(std::make_tuple(i, j, 1.0));
        }
    }
# endif
    return to_matrix(simil_tri, I, J, false, true);
}

/***R
system.time(
    out <- cpp_mask(sample(1:100, 10000, replace = TRUE), sample(1:100, 10000, replace = TRUE))
)
*/
