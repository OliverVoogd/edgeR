#include "glm.h"
#include "objects.h"

/* Different initialization methods for the Levenberg coefficients */

const char side='L';
const char trans_ormqr='T';
const char uplo='U';
const char trans_trtrs='N';
const char diag='N';
const int unity=1;

struct QRdecomposition {
    QRdecomposition(int nrows, int ncoefs, const double* curX) : NR(nrows), NC(ncoefs),
            X(curX), Xcopy(NR*NC), tau(NC), effects(NR), weights(NR), 
            lwork_geqrf(-1), lwork_ormqr(-1) {

        // Setting up the workspace for dgeqrf.
        double tmpwork;
        F77_CALL(dgeqrf)(&NR, &NC, Xcopy.data(), &NR, tau.data(), &tmpwork, &lwork_geqrf, &info);

        // Loading up the optimal WORK.
        lwork_geqrf=tmpwork+0.5;
        if (lwork_geqrf < 1) { lwork_geqrf = 1; }
        work_geqrf.resize(lwork_geqrf);

        // Repeating for dormqr
        F77_CALL(dormqr)(&side, &trans_ormqr, &NR, &unity, &NC, Xcopy.data(), &NR, tau.data(), 
                effects.data(), &NR, &tmpwork, &lwork_ormqr, &info);
        lwork_ormqr=tmpwork+0.5;
        if (lwork_ormqr < 1) { lwork_ormqr = 1; }
        work_ormqr.resize(lwork_ormqr);
        return;
    }

    void store_weights(const double* w) {
        if (w==NULL) {
            std::fill(weights.begin(), weights.end(), 1);
        } else {
            for (int row=0; row<NR; ++row) {
                weights[row]=std::sqrt(w[row]);
            }
        }
        return;
    }

    void decompose() {
        auto xcIt=Xcopy.begin();
        std::copy(X, X + Xcopy.size(), xcIt);
        for (int coef=0; coef<NC; ++coef) {
            for (int lib=0; lib<NR; ++lib) {
                (*xcIt)*=weights[lib];
                ++xcIt;
            }
        }
 
        F77_CALL(dgeqrf)(&NR, &NC, Xcopy.data(), &NR, tau.data(), work_geqrf.data(), &lwork_geqrf, &info);
        if (info) {
            throw std::runtime_error("QR decomposition failed");
        }
       return;
    } 

    void solve(const double * y) {
        for (int row=0; row<NR; ++row) {
            effects[row]=y[row]*weights[row];
        }

        F77_CALL(dormqr)(&side, &trans_ormqr, &NR, &unity, &NC, Xcopy.data(), &NR, tau.data(), 
                effects.data(), &NR, work_ormqr.data(), &lwork_ormqr, &info);
        if (info) {
            throw std::runtime_error("Q**T multiplication failed");
        }

        F77_CALL(dtrtrs)(&uplo, &trans_trtrs, &diag, &NC, &unity, Xcopy.data(), &NR, effects.data(), &NR, &info);
        if (info) {
            throw std::runtime_error("failed to solve the triangular system");
        }

        return;
    }

    const int NR, NC;
    const double* X;
    std::vector<double> Xcopy, tau, effects, weights, work_geqrf, work_ormqr;
    int lwork_geqrf, lwork_ormqr, info;
};

SEXP get_levenberg_start(SEXP y, SEXP offset, SEXP disp, SEXP weights, SEXP design, SEXP use_null) {
    BEGIN_RCPP
    any_numeric_matrix counts(y);
    const int num_tags=counts.get_nrow();
    const int num_libs=counts.get_ncol();

    Rcpp::NumericMatrix X=check_design_matrix(design, num_libs);
    const int num_coefs=X.ncol();
    QRdecomposition QR(num_libs, num_coefs, X.begin());

    // Initializing pointers to the assorted features.
    compressed_matrix allo=check_CM_dims(offset, num_tags, num_libs, "offset", "count");
    compressed_matrix alld=check_CM_dims(disp, num_tags, num_libs, "dispersion", "count");
    compressed_matrix allw=check_CM_dims(weights, num_tags, num_libs, "weight", "count");

    // Determining what type of algorithm is to be used.
    bool null_method=check_logical_scalar(use_null, "'use_null' specification");

    Rcpp::NumericMatrix output(num_tags, num_coefs);
    std::vector<double> current(num_libs);
    if (null_method) {
        QR.store_weights(NULL);
        QR.decompose();

        for (int tag=0; tag<num_tags; ++tag) {
            counts.fill_row(tag, current.data());
            const double* dptr=alld.get_row(tag);
            const double* optr=allo.get_row(tag);
            const double* wptr=allw.get_row(tag);
               
            // Computing weighted average of the count:library size ratios.
            double sum_weight=0, sum_exprs=0;
            for (int lib=0; lib<num_libs; ++lib) {
                const double curN=std::exp(optr[lib]);
                const double curweight=wptr[lib]*curN/(1 + dptr[lib] * curN);
                sum_exprs += current[lib] * curweight / curN;
                sum_weight += curweight;
            }
            std::fill(current.begin(), current.end(), std::log(sum_exprs/sum_weight));

            // Performing the QR decomposition and taking the solution.
            QR.solve(current.data());
            auto curout=output.row(tag);
            std::copy(QR.effects.begin(), QR.effects.begin()+num_coefs, curout.begin());
        }
    } else {
        const bool weights_are_the_same=allw.is_row_repeated();
        if (weights_are_the_same && num_tags) {
            QR.store_weights(allw.get_row(0));
            QR.decompose();
        }

        // Finding the delta.
        double delta=0;
        if (counts.is_data_integer()) {
            Rcpp::IntegerMatrix imat=counts.get_raw_int();
            delta=*std::max_element(imat.begin(), imat.end());
        } else {
            Rcpp::NumericMatrix dmat=counts.get_raw_dbl();
            delta=*std::max_element(dmat.begin(), dmat.end());
        }
        delta=std::min(delta, 1.0/6);

        for (int tag=0; tag<num_tags; ++tag) {
            if (!weights_are_the_same) {
                QR.store_weights(allw.get_row(tag));
                QR.decompose();
            }
            counts.fill_row(tag, current.data());
            const double* optr=allo.get_row(tag);
          
            // Computing normalized log-expression values.
            for (int lib=0; lib<num_libs; ++lib) {
                current[lib]=std::log(std::max(delta, current[lib])) - optr[lib];
            }

            // Performing the QR decomposition and taking the solution.
            QR.solve(current.data());
            auto curout=output.row(tag);
            std::copy(QR.effects.begin(), QR.effects.begin()+num_coefs, curout.begin());
        }
    }

    return output;
    END_RCPP
}
