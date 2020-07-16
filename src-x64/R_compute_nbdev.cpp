#include "utils.h"
#include "glm.h"
#include "objects.h"

SEXP compute_nbdev (SEXP y, SEXP mu, SEXP phi, SEXP weights, SEXP dosum) {
    BEGIN_RCPP
    any_numeric_matrix counts(y);
    const int num_tags=counts.get_nrow();
    const int num_libs=counts.get_ncol();
    std::vector<double> current(num_libs);

    // Setting up means.
    Rcpp::NumericMatrix fitted(mu);
    if (fitted.nrow()!=num_tags || fitted.ncol()!=num_libs) {
        throw std::runtime_error("dimensions of count and fitted value matrices are not equal");
    }

    // Setting up dispersions.
    compressed_matrix alld=check_CM_dims(phi, num_tags, num_libs, "dispersion", "count");

    // Seeing if we have to sum things together.
    bool sumtogether=check_logical_scalar(dosum, "summation specifier");

    if (sumtogether) {
        // Setting up weights.
        compressed_matrix allw(weights);

        Rcpp::NumericVector output(num_tags);
        for (int tag=0; tag<num_tags; ++tag) {
            counts.fill_row(tag, current.data());
            const double* dptr=alld.get_row(tag);
            const double* wptr=allw.get_row(tag);

            Rcpp::NumericMatrix::Row curmeans=fitted.row(tag);
            double& current_sumdev=output[tag];
            auto cmIt=curmeans.begin(); 
            for (int lib=0; lib<num_libs; ++lib, ++cmIt) {
                current_sumdev += compute_unit_nb_deviance(current[lib], *cmIt, dptr[lib]) * wptr[lib];
            }
        }

        return output;
    } else {
        Rcpp::NumericMatrix output(num_tags, num_libs);
        for (int tag=0; tag<num_tags; ++tag) {
            counts.fill_row(tag, current.data());
            const double* dptr=alld.get_row(tag);

            auto curmeans=fitted.row(tag);
            auto cmIt=curmeans.begin();
            auto outvals=output.row(tag);
            auto ovIt=outvals.begin();

            for (int lib=0; lib<num_libs; ++lib, ++ovIt, ++cmIt) {
                (*ovIt) = compute_unit_nb_deviance(current[lib], *cmIt, dptr[lib]);
            }
       }

        return output;
    }
    END_RCPP
}
