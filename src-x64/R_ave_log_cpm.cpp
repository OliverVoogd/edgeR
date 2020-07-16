#include "glm.h"
#include "add_prior.h"
#include "objects.h"

SEXP ave_log_cpm(SEXP y, SEXP offset, SEXP prior, SEXP disp, SEXP weights, SEXP max_iterations, SEXP tolerance) {
    BEGIN_RCPP

    any_numeric_matrix counts(y);
    const int num_tags=counts.get_nrow();
    const int num_libs=counts.get_ncol();
    std::vector<double> current(num_libs);

    add_prior AP(prior, offset, true, true);
    check_AP_dims(AP, num_tags, num_libs, "count");
    compressed_matrix alld=check_CM_dims(disp, num_tags, num_libs, "dispersion", "count");
    compressed_matrix allw=check_CM_dims(weights, num_tags, num_libs, "weight", "count");
    
    // GLM fitting specifications
    int maxit=check_integer_scalar(max_iterations, "maximum iterations");
    double tol=check_numeric_scalar(tolerance, "tolerance");

    // Returning average log-cpm
    Rcpp::NumericVector output(num_tags);
    for (int tag=0; tag<num_tags; ++tag) {
        counts.fill_row(tag, current.data());
           
        // Adding the current set of priors.
        AP.compute(tag);
        const double* offptr=AP.get_offsets();
        const double* pptr=AP.get_priors();
        for (int lib=0; lib<num_libs; ++lib) {
            current[lib]+=pptr[lib];
        }

        // Fitting a one-way layout.
        const double* dptr=alld.get_row(tag);
        const double* wptr=allw.get_row(tag);
        auto fit=glm_one_group(num_libs, current.data(), offptr, dptr, wptr, maxit, tol, NA_REAL);
        output[tag]=(fit.first + LNmillion)/LNtwo;
    }
    
    return output;
    END_RCPP
}
