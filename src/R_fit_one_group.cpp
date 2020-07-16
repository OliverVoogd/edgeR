#include "glm.h"
#include "objects.h"

template <typename T>
bool is_array_equal_to(const T* x, const int n, const bool rep, const T& v) {
    if (rep) {
        return (n>0 && x[0]==v);
    } else {
        for (int i=0; i<n; ++i) {
            if (x[i]!=v) { return false; }
        }
        return true;
    }
}

SEXP fit_one_group (SEXP y, SEXP offsets, SEXP disp, SEXP weights, SEXP max_iterations, SEXP tolerance, SEXP beta) {
    BEGIN_RCPP
    any_numeric_matrix counts(y);
    const int num_tags=counts.get_nrow();
    const int num_libs=counts.get_ncol();
    std::vector<double> current(num_libs);

    // Setting up assorted input matrices.
    compressed_matrix allo=check_CM_dims(offsets, num_tags, num_libs, "offset", "count");
    compressed_matrix alld=check_CM_dims(disp, num_tags, num_libs, "dispersion", "count");
    compressed_matrix allw=check_CM_dims(weights, num_tags, num_libs, "weight", "count");

    // Setting up the beta object.
    Rcpp::NumericVector Beta(beta);
    if (Beta.size()!=num_tags) {
        throw std::runtime_error("length of beta vector should be equal to number of genes");
    }

    // Setting up scalars.
    int maxit=check_integer_scalar(max_iterations, "maximum iterations");
    double tol=check_numeric_scalar(tolerance, "tolerance");

    // Setting up beta for output.
    Rcpp::NumericVector out_beta(num_tags);
    Rcpp::LogicalVector out_conv(num_tags);

    // Preparing for possible Poisson sums.
    bool disp_is_zero=true, weight_is_one=true;
    double sum_lib=0;
    if (allo.is_row_repeated() && num_tags) {
        const double* optr=allo.get_row(0);
        for (int lib=0; lib<num_libs; ++lib) {
            sum_lib+=std::exp(optr[lib]);
        }
     }
    if (alld.is_row_repeated() && num_tags) {
        const double* dptr=alld.get_row(0);
        disp_is_zero=is_array_equal_to<double>(dptr, num_libs, alld.is_col_repeated(), 0);
    }
    if (allw.is_row_repeated() && num_tags) {
        const double* wptr=allw.get_row(0);
        weight_is_one=is_array_equal_to<double>(wptr, num_libs, allw.is_col_repeated(), 1);
    }

    // Iterating through tags and fitting.
	for (int tag=0; tag<num_tags; ++tag) {
        counts.fill_row(tag, current.data());
        const double* optr=allo.get_row(tag);
        const double* wptr=allw.get_row(tag);
        const double* dptr=alld.get_row(tag);

        // Checking for the Poisson special case with all-unity weights and all-zero dispersions.
        if (!alld.is_row_repeated()) {
            disp_is_zero=is_array_equal_to<double>(dptr, num_libs, alld.is_col_repeated(), 0);
        }
        if (!allw.is_row_repeated()) {
            weight_is_one=is_array_equal_to<double>(wptr, num_libs, allw.is_col_repeated(), 1);
        }

        if (disp_is_zero && weight_is_one) {
            if (!allo.is_row_repeated()) {
                // Only recalculate sum of library sizes if it has changed.
                sum_lib=0;
                for (int lib=0; lib<num_libs; ++lib) { sum_lib+=std::exp(optr[lib]); }
            }

            double sum_counts=std::accumulate(current.begin(), current.end(), 0.0);
            if (sum_counts==0) {
                out_beta[tag]=R_NegInf;
            } else {
                out_beta[tag]=std::log(sum_counts/sum_lib);
            }
            out_conv[tag]=true;
        } else {
            // Otherwise going through NR iterations.
            std::pair<double, bool> out=glm_one_group(num_libs, current.data(), optr, dptr, wptr, maxit, tol, Beta[tag]);
            out_beta[tag]=out.first;
            out_conv[tag]=out.second;
        }
	}

	// Returning everything as a list.
    return Rcpp::List::create(out_beta, out_conv);
    END_RCPP
}
