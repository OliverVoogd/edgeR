#include "objects.h"

/* Checks whether the variance is below the Poisson bound. */

SEXP check_poisson_bound (SEXP fitted, SEXP disp, SEXP s2) {
    BEGIN_RCPP
    Rcpp::NumericMatrix Fitted(fitted);
    const int num_tags=Fitted.nrow();
    const int num_libs=Fitted.ncol();

    compressed_matrix alld=check_CM_dims(disp, num_tags, num_libs, "NB dispersion", "fitted value");
    compressed_matrix alls=check_CM_dims(s2, num_tags, num_libs, "QL dispersion", "fitted value");

    Rcpp::LogicalVector output(num_tags);
    for (int tag=0; tag<num_tags; ++tag) {
        int& below_bound=output[tag];
        const double* dptr=alld.get_row(tag);
        const double* sptr=alls.get_row(tag);

        auto current=Fitted.row(tag);
        for (const auto& val : current) { 
            if ((val * (*dptr) + 1) * (*sptr) < 1) {
                below_bound=1;
                break;
            }
            ++dptr;
            ++sptr;
        }
    }
    
    return output;
    END_RCPP
}
