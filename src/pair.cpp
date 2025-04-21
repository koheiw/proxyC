#include "proxyc.h"
#include "dev.h"
using namespace proxyc;
using namespace arma;

double simil_cosine(colvec& col_i, colvec& col_j) {
    return accu(sum(col_i % col_j) / sqrt(sum(square(col_i)) * sum(square(col_j))));
}

double simil_correlation(colvec& col_i, colvec& col_j) {
    double sd_i = stddev(col_i, 1);
    double sd_j = stddev(col_j, 1);
    double v1 = accu(col_i.t() * col_j);
    double v2 = mean(col_i) * mean(col_j) * col_i.n_rows;
    // Rcout << "v1 - v2 = " << v1 - v2 << "\n";
    // Rcout << "col_i.n_rows = " << col_i.n_rows << "\n";
    // Rcout << "sd_i * sd_j = " << sd_i * sd_j << "\n";
    // Rcout << "result = " << ((v1 - v2) / col_i.n_rows) / (sd_i * sd_j) << "\n";
    return ((v1 - v2) / col_i.n_rows) / (sd_i * sd_j);

}

double simil_ejaccard(colvec& col_i, colvec& col_j, double weight = 1) {
    double e = accu(col_i % col_j);
    if (e == 0) return 0;
    return e / (accu(pow(col_i, weight)) + accu(pow(col_j, weight)) - e);
}

double simil_fjaccard(colvec& col_i, colvec& col_j) {
    if (any(col_i < 0) || any(1.0 < col_i) ||
        any(col_j < 0) || any(1.0 < col_j))
        return std::numeric_limits<double>::quiet_NaN();
    ucolvec l = (col_i <= col_j);
    colvec min = (col_i % l) + (col_j % (1 - l));
    colvec max = (col_i % (1 - l)) + (col_j % l);
    // Rcout << "min\n" << min.t();
    // Rcout << "max\n" << max.t();
    return accu(min) / accu(max);
}

double simil_dice(colvec& col_i, colvec& col_j) {
    double e = accu(col_i % col_j);
    if (e == 0) return 0;
    return (2 * e) / (accu(col_i) + accu(col_j));
}

double simil_edice(colvec& col_i, colvec& col_j) {
    double e = accu(col_i % col_j);
    if (e == 0) return 0;
    return (2 * e) / (accu(square(col_i)) + accu(square(col_j)));
}

double simil_hamann(colvec& col_i, colvec& col_j) {
    double e = accu(col_i == col_j);
    double u = col_i.n_rows - e;
    return (e - u) / (e + u);
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
        return std::numeric_limits<double>::quiet_NaN();
    mat m = join_rows(col_i, col_j) + smooth;
    m = m / accu(m);
    mat e = sum(m, 1) * sum(m, 0);
    mat d = pow(m - e, 2) / e;
    return accu(d);
}

double dist_kullback(colvec& col_i, colvec& col_j, double smooth) {
    if (smooth == 0 && (any(col_i == 0) || any(col_j == 0)))
        return std::numeric_limits<double>::quiet_NaN();
    double s1 = accu(col_i) + smooth * col_i.n_rows;
    double s2 = accu(col_j) + smooth * col_j.n_rows;
    colvec p1 = (col_i + smooth) / s1;
    colvec p2 = (col_j + smooth) / s2;
    return accu(trans(p2) * log(p2 / p1));
}

double dist_jeffreys(colvec& col_i, colvec& col_j, double smooth) {
    if (smooth == 0 && (any(col_i == 0) || any(col_j == 0)))
        return std::numeric_limits<double>::quiet_NaN();
    double s1 = accu(col_i) + smooth * col_i.n_rows;
    double s2 = accu(col_j) + smooth * col_j.n_rows;
    colvec p1 = (col_i + smooth) / s1;
    colvec p2 = (col_j + smooth) / s2;
    return accu(trans(p2 - p1) * log(p2 / p1));
}

double dist_jensen(colvec& col_i, colvec& col_j, double smooth) {
    if (smooth == 0 && (any(col_i == 0) || any(col_j == 0)))
        return std::numeric_limits<double>::quiet_NaN();
    double s1 = accu(col_i) + smooth * col_i.n_rows;
    double s2 = accu(col_j) + smooth * col_j.n_rows;
    colvec p1 = (col_i + smooth) / s1;
    colvec p2 = (col_j + smooth) / s2;
    colvec m = (p1 + p2) / 2;
    return (accu(trans(p1) * log(p1 / m)) + accu(trans(p2) * log(p2 / m))) / 2;
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

void proxy_pair(const uword i,
                const sp_mat& mt1, const sp_mat& mt2, const sp_mat& mask,
                Triplets& simil_tri,
                const int method,
                const unsigned int rank, const double limit, const bool symm,
                const bool diag, const double weight, const double smooth,
                const bool drop0, const bool use_nan, const bool use_mask,
                const int digits) {

    uword nrow = mt1.n_rows;
    uword ncol = mt1.n_cols;

    colvec col_i(nrow);
    colvec col_j(nrow);
    double simil = 0;
    std::vector<double> simils;
    col_i = mt2.col(i);
    if (diag) {
        simils.reserve(1);
    } else {
        simils.reserve(ncol);
    }
    colvec mask_i;
    if (use_mask)
        mask_i = colvec(mask.col(i));
    for (uword j = 0; j < ncol; j++) {
        if (use_mask && mask_i.at(j) == 0) {
            simils.push_back(std::numeric_limits<double>::quiet_NaN());
            continue;
        };
        if (diag && j != i) continue;
        if (symm && j > i) continue;
        col_j = mt1.col(j);
        switch (method) {
        case 1:
            simil = simil_cosine(col_i, col_j);
            break;
        case 2:
            simil = simil_correlation(col_i, col_j);
            simil = replace_inf(simil);
            break;
        case 3:
            simil = dist_euclidean(col_i, col_j);
            break;
        case 4:
            simil = simil_dice(col_i, col_j);
            break;
        case 5:
            simil = simil_edice(col_i, col_j);
            break;
        case 6:
            simil = simil_hamann(col_i, col_j);
            break;
        case 7:
            simil = simil_matching(col_i, col_j);
            break;
        case 8:
            simil = simil_faith(col_i, col_j);
            break;
        case 9:
            simil = simil_ejaccard(col_i, col_j, weight);
            break;
        case 10:
            simil = simil_fjaccard(col_i, col_j);
            break;
        case 11:
            simil = dist_chisquare(col_i, col_j, smooth);
            break;
        case 12:
            simil = dist_kullback(col_i, col_j, smooth);
            break;
        case 13:
            simil = dist_manhattan(col_i, col_j);
            break;
        case 14:
            simil = dist_maximum(col_i, col_j);
            break;
        case 15:
            simil = dist_canberra(col_i, col_j);
            break;
        case 16:
            simil = dist_minkowski(col_i, col_j, weight);
            break;
        case 17:
            simil = dist_hamming(col_i, col_j);
            break;
        case 18:
            simil = dist_jeffreys(col_i, col_j, smooth);
            break;
        case 19:
            simil = dist_jensen(col_i, col_j, smooth);
            break;
        }
        //Rcout << "simil=" << simil << "\n";
        simils.push_back(simil);
    }
    simils = round(simils, digits);
    double l = get_limit(simils, rank, limit);
    for (std::size_t k = 0; k < simils.size(); k++) {
        double s = simils[k];
        if (drop0 && s == 0) continue;
        if (s >= l || (use_nan && std::isnan(s))) {
            if (diag) {
                simil_tri.push_back(std::make_tuple(i, i, s));
            } else {
                simil_tri.push_back(std::make_tuple(k, i, s));
            }
        }
    }
}

// [[Rcpp::export]]
S4 cpp_pair(arma::sp_mat& mt1,
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
    out <- cpp_pair(mt1, mt2, 1, mask, 100, symm = TRUE, use_mask = TRUE, drop0 = TRUE)
)

*/
