#include "proxyc.h"
//#include "dev.h"
using namespace proxyc;
using namespace arma;

namespace proxyc{

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

    inline rowvec stddev(const sp_mat& mt, const int norm_type) {
        rowvec v(mt.n_cols, fill::zeros);
        if (mt.is_empty()) return(v);
        for (uword i = 0; i < mt.n_cols; i++) {
            v[i] = arma::stddev(colvec(mt.col(i)), norm_type);
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

    inline rowvec mean(const sp_mat& mt) {
        rowvec v(mt.n_cols, fill::zeros);
        if (mt.is_empty()) return(v);
        for (uword i = 0; i < mt.n_cols; i++) {
            v[i] = arma::mean(colvec(mt.col(i)));
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

    inline rowvec max(const sp_mat& mt) {
        rowvec v(mt.n_cols, fill::zeros);
        if (mt.is_empty()) return(v);
        for (uword i = 0; i < mt.n_cols; i++) {
            v[i] = arma::max(colvec(mt.col(i)));
        }
        return(v);
    }

    inline rowvec max(const mat& mt) {
        rowvec v(mt.n_cols, fill::zeros);
        if (mt.is_empty()) return(v);
        for (uword i = 0; i < mt.n_cols; i++) {
            v[i] = arma::max(colvec(mt.col(i)));
        }
        return(v);
    }

    inline rowvec min(const sp_mat& mt) {
        rowvec v(mt.n_cols, fill::zeros);
        if (mt.is_empty()) return(v);
        for (uword i = 0; i < mt.n_cols; i++) {
            v[i] = arma::min(colvec(mt.col(i)));
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

    inline urowvec maxrow(const sp_mat& mt) {
        urowvec v(mt.n_cols, fill::zeros);
        if (mt.is_empty()) return(v);
        for (uword i = 0; i < mt.n_cols; i++) {
            v[i] = mt.col(i).index_max();
        }
        return(v);
    }

    inline urowvec maxrow(const mat& mt) {
        urowvec v(mt.n_cols, fill::zeros);
        if (mt.is_empty()) return(v);
        for (uword i = 0; i < mt.n_cols; i++) {
            v[i] = mt.col(i).index_max();
        }
        return(v);
    }

}
