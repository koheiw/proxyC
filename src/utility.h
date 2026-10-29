#include "proxyc.h"
//#include "dev.h"
using namespace proxyc;
using namespace arma;

namespace proxyc{

    // nnz
    inline urowvec nnz(const sp_mat& mt) {
        urowvec v(mt.n_cols, fill::zeros);
        if (mt.is_empty()) return v;
        for (uword i = 0; i < mt.n_cols; i++) {
            v[i] = arma::sum(colvec(mt.col(i)) != 0);
        }
        return v;
    }

    inline urowvec nnz(const mat& mt) {
        urowvec v(mt.n_cols, fill::zeros);
        if (mt.is_empty()) return v;
        for (uword i = 0; i < mt.n_cols; i++) {
            v[i] = arma::accu(mt.col(i) != 0);
        }
        return v;
    }

    // stddev
    inline rowvec stddev(const sp_mat& mt, const int norm_type,
                         const bool use_nan = false) {
        rowvec v(mt.n_cols, fill::zeros);
        if (mt.is_empty()) return v;
        colvec e(mt.n_rows);
        for (uword i = 0; i < mt.n_cols; i++) {
            e = colvec(mt.col(i));
            if (use_nan && !arma::any(e)) {
                v[i] = arma::datum::nan;
            } else {
                v[i] = arma::stddev(e, norm_type);
            }
        }
        return v;
    }

    inline rowvec stddev(const mat& mt, const int norm_type) {
        return(arma::stddev(mt, norm_type, 0));
    }

    // mean
    inline rowvec mean(const sp_mat& mt, const bool use_nan = false) {
        rowvec v(mt.n_cols, fill::zeros);
        if (mt.is_empty()) return v;
        colvec e(mt.n_rows);
        for (uword i = 0; i < mt.n_cols; i++) {
            e = colvec(mt.col(i));
            if (use_nan && !arma::any(e)) {
                v[i] = arma::datum::nan;
            } else {
                v[i] = arma::mean(e);
            }
        }
        return v;
    }

    inline rowvec mean(const mat& mt) {
        return(arma::mean(mt, 0));
    }

    // max
    inline rowvec max(const sp_mat& mt, const bool use_nan = false) {
        rowvec v(mt.n_cols, fill::zeros);
        if (mt.is_empty()) return v;
        colvec e(mt.n_rows);
        for (uword i = 0; i < mt.n_cols; i++) {
            e = colvec(mt.col(i));
            if (use_nan && !arma::any(e)) {
                v[i] = arma::datum::nan;
            } else {
                v[i] = arma::max(e);
            }
        }
        return v;
    }

    inline rowvec max(const mat& mt) {
        return(arma::max(mt, 0));
    }

    // min
    inline rowvec min(const sp_mat& mt, const bool use_nan = false) {
        rowvec v(mt.n_cols, fill::zeros);
        if (mt.is_empty()) return v;
        colvec e(mt.n_rows);
        for (uword i = 0; i < mt.n_cols; i++) {
            e = colvec(mt.col(i));
            if (use_nan && !arma::any(e)) {
                v[i] = arma::datum::nan;
            } else {
                v[i] = arma::min(e);
            }
        }
        return v;
    }

    inline rowvec min(const mat& mt) {
        return(arma::min(mt, 0));
    }

    // maxrow
    inline rowvec maxrow(const sp_mat& mt, const int mode = 1,
                         const bool use_nan = false) {
        rowvec v(mt.n_cols, fill::zeros);
        if (mt.is_empty()) return v;
        colvec e(mt.n_rows);
        for (uword i = 0; i < mt.n_cols; i++) {
            e = colvec(mt.col(i));
            if (use_nan && !arma::any(e)) {
                v[i] = arma::datum::nan;
            } else {
                uvec idx;
                switch (mode) {
                case 1: // random
                    idx = shuffle(find(e == max(e)));
                    break;
                case 2: // first
                    idx = find(e == max(e), 1);
                    break;
                case 3: // last
                    idx = find(e == max(e), 1, "last");
                    break;
                }
                if (idx.size() == 0) {
                    v[i] = arma::datum::nan;
                } else {
                    v[i] = idx(0);
                }
            }
        }
        return v;
    }

    inline rowvec maxrow(const mat& mt, const int mode = 1) {
        rowvec v(mt.n_cols, fill::zeros);
        if (mt.is_empty()) return v;
        colvec e(mt.n_rows);
        for (uword i = 0; i < mt.n_cols; i++) {
            e = colvec(mt.col(i));
            uvec idx;
            switch (mode) {
            case 1: // random
                idx = shuffle(find(e == max(e)));
                break;
            case 2: // first
                idx = find(e == max(e), 1);
                break;
            case 3: // last
                idx = find(e == max(e), 1, "last");
                break;
            }
            if (idx.size() == 0) {
                v[i] = arma::datum::nan;
            } else {
                v[i] = idx(0);
            }
        }
        return v;
    }

}
