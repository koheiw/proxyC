#define ARMA_NO_DEBUG
#include "armadillo.h"
#include "proxyc.h"
#include "dev.h"
using namespace proxyc;
using namespace arma;

double simil_cosine(colvec& col_i, colvec& col_j) {
    return accu(sum(col_i % col_j) / sqrt(sum(square(col_i)) * sum(square(col_j))));
}

double simil_correlation(colvec& col_i, colvec& col_j) {
    double v1 = accu(col_i.t() * col_j);
    double v2 = mean(col_i) * mean(col_j) * col_i.n_rows;
    return ((v1 - v2) / col_i.n_rows) / (stddev(col_i, 1) * stddev(col_j, 1));
}

double simil_ejaccard(colvec& col_i, colvec& col_j, double weight = 1) {
    double e = accu(col_i % col_j);
    return e / (accu(pow(col_i, weight)) + accu(pow(col_j, weight)) - e);
}

double simil_edice(colvec& col_i, colvec& col_j, double weight = 1) {
    double e = accu(col_i % col_j);
    return (2 * e) / (accu(pow(col_i, weight)) + accu(pow(col_j, weight)));
}

double simil_hamman(colvec& col_i, colvec& col_j, double weight = 1) {
    double e = accu(col_i == col_j);
    double u = col_i.n_rows - e;
    return (e - (u * weight)) / (e + u);
}

double simil_faith(colvec& col_i, colvec& col_j) {
    double t = accu(col_i == 1 && col_j == 1);
    double f = accu(col_i == 0 && col_j == 0);
    double n = col_i.n_rows;
    return (t + (f / 2)) / n;
}

double simil_matching(colvec& col_i, colvec& col_j) {
    uvec m = col_i == col_j;
    double n = m.n_rows;
    return accu(m) / n;
}

double dist_euclidean(colvec& col_i, colvec& col_j) {
    return sqrt(accu(square(col_i - col_j)));
}

double dist_chisquare(colvec& col_i, colvec& col_j, double smooth) {
    if (smooth == 0 && (any(col_i == 0) || any(col_j == 0)))
        return 0;
    mat m = join_rows(col_i, col_j) + smooth;
    m = m / accu(m);
    mat e = sum(m, 1) * sum(m, 0);
    mat d = pow(m - e, 2) / e;
    return accu(d);
}

double dist_kullback(colvec& col_i, colvec& col_j, double smooth) {
    if (smooth == 0 && (any(col_i == 0) || any(col_j == 0)))
        return 0;
    double s1 = accu(col_i) + smooth * col_i.n_rows;
    double s2 = accu(col_j) + smooth * col_j.n_rows;
    colvec p1 = (col_i + smooth) / s1;
    colvec p2 = (col_j + smooth) / s2;
    return accu(trans(p2) * log(p2 / p1));
}

double dist_manhattan(colvec& col_i, colvec& col_j) {
    return accu(abs(col_i - col_j));
}

double dist_hamming(colvec& col_i, colvec& col_j) {
    return accu(col_i != col_j);
}

double dist_maximum(colvec& col_i, colvec& col_j) {
    return accu(max(abs(col_i - col_j)));
}

double dist_canberra(colvec& col_i, colvec& col_j) {
    colvec b = abs(col_i - col_j);
    colvec m = abs(col_i) + abs(col_j);
    uvec nz = find(m > 0);
    return accu(b(nz) / m(nz));
}

double dist_minkowski(colvec& col_i, colvec& col_j, double p = 1) {
    return pow(accu(pow(abs(col_i - col_j), p)), 1 / p);
}

struct pairWorker : public Worker {

    const sp_mat& mt1; // input
    const sp_mat& mt2; // input
    Triplets& simil_tri; // output
    const int method;
    const unsigned int rank;
    const double limit;
    const bool symm;
    const bool diag;
    const double weight;
    const double smooth;
    const bool drop0;

    pairWorker(const sp_mat& mt1_, const sp_mat& mt2_, Triplets& simil_tri_,
               const int method_,
               const unsigned int rank_, const double limit_, const bool symm_,
               const bool diag_, const double weight_, const double smooth_, const bool drop0_) :
               mt1(mt1_), mt2(mt2_), simil_tri(simil_tri_),
               method(method_), rank(rank_), limit(limit_), symm(symm_), diag(diag_),
               weight(weight_), smooth(smooth_), drop0(drop0_) {}

    void operator()(std::size_t begin, std::size_t end) {

        arma::uword nrow = mt1.n_rows;
        arma::uword ncol = mt1.n_cols;

        double simil = 0;
        std::vector<double> simils;
        colvec col_i(nrow);
        colvec col_j(nrow);
        for (uword i = begin; i < end; i++) {
            col_i = mt2.col(i);
            if (diag) {
                simils.reserve(1);
            } else {
                simils.reserve(ncol);
            }
            for (uword j = 0; j < ncol; j++) {
                if (diag && j != i) continue;
                if (symm && j > i) continue;
                col_j = mt1.col(j);
                switch (method) {
                case 1:
                    simil = simil_cosine(col_i, col_j);
                    break;
                case 2:
                    simil = simil_correlation(col_i, col_j);
                    break;
                case 3:
                    simil = simil_ejaccard(col_i, col_j, weight);
                    break;
                case 4:
                    simil = simil_edice(col_i, col_j, weight);
                    break;
                case 5:
                    simil = simil_hamman(col_i, col_j, weight);
                    break;
                case 6:
                    simil = simil_matching(col_i, col_j);
                    break;
                case 7:
                    simil = simil_faith(col_i, col_j);
                    break;
                case 8:
                    simil = dist_euclidean(col_i, col_j);
                    break;
                case 9:
                    simil = dist_chisquare(col_i, col_j, smooth);
                    break;
                case 10:
                    simil = dist_kullback(col_i, col_j, smooth);
                    break;
                case 11:
                    simil = dist_manhattan(col_i, col_j);
                    break;
                case 12:
                    simil = dist_maximum(col_i, col_j);
                    break;
                case 13:
                    simil = dist_canberra(col_i, col_j);
                    break;
                case 14:
                    simil = dist_minkowski(col_i, col_j, weight);
                    break;
                case 15:
                    simil = dist_hamming(col_i, col_j);
                    break;
                }
                //Rcout << "simil=" << simil << "\n";
                simils.push_back(simil);
            }
            double l = get_limit(simils, rank, limit);
            for (std::size_t k = 0; k < simils.size(); k++) {
                if (simils[k] >= l) {
                    if (drop0 && simils[k] == 0) continue;
                    if (diag) {
                        simil_tri.push_back(std::make_tuple(i, i, simils[k]));
                    } else {
                        simil_tri.push_back(std::make_tuple(k, i, simils[k]));
                    }
                }
            }
            simils.clear();
        }
    }
};

// [[Rcpp::export]]
S4 cpp_pair(arma::sp_mat& mt1,
            arma::sp_mat& mt2,
            const int method,
            unsigned int rank,
            double limit = -1.0,
            double weight = 1.0,
            double smooth = 0,
            bool symm = false,
            bool diag = false,
            bool drop0 = false) {

    if (mt1.n_rows != mt2.n_rows)
        throw std::range_error("Invalid matrix objects");

    uword ncol1 = mt1.n_cols;
    uword ncol2 = mt2.n_cols;
    if (rank < 1) rank = 1;
    symm = symm && method != 10 && rank == ncol2; // exception for kullback

    //dev::Timer timer;
    //dev::start_timer("Compute similarity", timer);
    Triplets simil_tri;
    pairWorker proxy_pair(mt1, mt2, simil_tri, method, rank, limit, symm, diag,
                          weight, smooth, drop0);
    parallelFor(0, ncol2, proxy_pair);
    //dev::stop_timer("Compute similarity", timer);

    return to_matrix(simil_tri, ncol1, ncol2, symm);

}

/***R
mt <- Matrix::rsparsematrix(1000, 500, 0.01)
microbenchmark::microbenchmark(
out <- cpp_pair(mt, mt, 6, 100, symm = TRUE, smooth = 0.0)
)
*/
