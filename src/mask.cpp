#include "proxyc.h"
#include "dev.h"
using namespace proxyc;
using namespace arma;

// [[Rcpp::export]]
S4 cpp_mask(IntegerVector v1, IntegerVector v2) {

    Triplets simil_tri;
    for (int i = 0; i < v1.size(); i++) {
        for (int j = 0; j < v2.size(); j++) {
            if (v1[i] == v2[j])
                simil_tri.push_back(std::make_tuple(i, j, 1.0));
        }
    }
    return to_matrix(simil_tri, v1.size(), v2.size(), false, true);

}

/***R
system.time(
    out <- cpp_mask(sample(1:100, 10000, replace = TRUE), sample(1:100, 10000, replace = TRUE))
)
*/
