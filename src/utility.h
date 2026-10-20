#include "proxyc.h"
//#include "dev.h"
using namespace proxyc;
using namespace arma;

namespace proxyc{

    // nnz
    inline rowvec nnz(const sp_mat& mt) {
        rowvec v(mt.n_cols, fill::zeros);
        if (mt.is_empty()) return(v);
        for (uword i = 0; i < mt.n_cols; i++) {
            v[i] = arma::sum(colvec(mt.col(i)) != 0);
        }
        return(v);
    }

    inline rowvec nnz(const mat& mt) {
        rowvec v(mt.n_cols, fill::zeros);
        if (mt.is_empty()) return(v);
        for (uword i = 0; i < mt.n_cols; i++) {
            v[i] = arma::sum(colvec(mt.col(i)) != 0);
        }
        return(v);
    }

    // stddev
    inline rowvec stddev(const sp_mat& mt, const int norm_type,
                         const bool use_nan = false) {
        rowvec v(mt.n_cols, fill::zeros);
        if (mt.is_empty()) return(v);
        colvec e(mt.n_rows);
        for (uword i = 0; i < mt.n_cols; i++) {
            e = colvec(mt.col(i));
            if (arma::any(e)) {
                v[i] = arma::stddev(e, norm_type);
            } else if (use_nan) {
                v[i] = arma::datum::nan;
            }
        }
        return(v);
    }

    inline rowvec stddev(const mat& mt, const int norm_type) {
        rowvec v(mt.n_cols, fill::zeros);
        if (mt.is_empty()) return(v);
        for (uword i = 0; i < mt.n_cols; i++) {
            v[i] = arma::stddev(colvec(mt.col(i)), norm_type);
        }
        return(v);
    }

    // mean
    inline rowvec mean(const sp_mat& mt, const bool use_nan = false) {
        rowvec v(mt.n_cols, fill::zeros);
        if (mt.is_empty()) return(v);
        colvec e(mt.n_rows);
        for (uword i = 0; i < mt.n_cols; i++) {
            e = colvec(mt.col(i));
            if (arma::any(e)) {
                v[i] = arma::mean(e);
            } else if (use_nan) {
                v[i] = arma::datum::nan;
            }
        }
        return(v);
    }

    inline rowvec mean(const mat& mt) {
        rowvec v(mt.n_cols, fill::zeros);
        if (mt.is_empty()) return(v);
        for (uword i = 0; i < mt.n_cols; i++) {
            v[i] = arma::mean(colvec(mt.col(i)));
        }
        return(v);
    }

    // max
    inline rowvec max(const sp_mat& mt, const bool use_nan = false) {
        rowvec v(mt.n_cols, fill::zeros);
        if (mt.is_empty()) return(v);
        colvec e(mt.n_rows);
        for (uword i = 0; i < mt.n_cols; i++) {
            e = colvec(mt.col(i));
            if (arma::any(e)) {
                v[i] = arma::max(e);
            } else if (use_nan) {
                v[i] = arma::datum::nan;
            }
        }
        return(v);
    }

    inline rowvec max(const mat& mt) {
        rowvec v(mt.n_cols, fill::nan);
        if (mt.is_empty()) return(v);
        for (uword i = 0; i < mt.n_cols; i++) {
            v[i] = arma::max(colvec(mt.col(i)));
        }
        return(v);
    }

    // min
    inline rowvec min(const sp_mat& mt, const bool use_nan = false) {
        rowvec v(mt.n_cols, fill::zeros);
        if (mt.is_empty()) return(v);
        colvec e(mt.n_rows);
        for (uword i = 0; i < mt.n_cols; i++) {
            e = colvec(mt.col(i));
            if (arma::any(e)) {
                v[i] = arma::min(e);
            } else if (use_nan) {
                v[i] = arma::datum::nan;
            }
        }
        return(v);
    }

    inline rowvec min(const mat& mt) {
        rowvec v(mt.n_cols, fill::zeros);
        if (mt.is_empty()) return(v);
        for (uword i = 0; i < mt.n_cols; i++) {
            v[i] = arma::min(colvec(mt.col(i)));
        }
        return(v);
    }

    // maxrow
    inline rowvec maxrow(const sp_mat& mt, const int mode = 1,
                         const bool use_nan = false) {
        rowvec v(mt.n_cols, fill::zeros);
        if (mt.is_empty()) return(v);
        colvec e(mt.n_rows);
        for (uword i = 0; i < mt.n_cols; i++) {
            e = colvec(mt.col(i));
            if (arma::any(e)) {
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
            } else if (use_nan) {
                v[i] = arma::datum::nan;
            }
        }
        return(v);
    }

    inline rowvec maxrow(const mat& mt, const int mode = 1) {
        rowvec v(mt.n_cols, fill::zeros);
        if (mt.is_empty()) return(v);
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
        return(v);
    }

}
