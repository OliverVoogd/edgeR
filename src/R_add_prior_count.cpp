#include "utils.h"
#include "add_prior.h"

/**** Adding a prior count to each observation. *******/

SEXP add_prior_count(SEXP y, SEXP offset, SEXP prior) {
    BEGIN_RCPP

    any_numeric_matrix input(y);
    const int num_tags=input.get_nrow();
    const int num_libs=input.get_ncol();

    Rcpp::NumericMatrix outmat(num_tags, num_libs);
    if (input.is_data_integer()) {
        auto tmp=input.get_raw_int();
        std::copy(tmp.begin(), tmp.end(), outmat.begin());
    } else {
        auto tmp=input.get_raw_dbl();
        std::copy(tmp.begin(), tmp.end(), outmat.begin());
    }

    add_prior AP(prior, offset, true, true);
    check_AP_dims(AP, num_tags, num_libs, "count");

    // Computing the adjusted library sizes, either as a vector or as a matrix.
    const bool same_prior=AP.same_across_rows();
    double* libptr=NULL;
    Rcpp::List output(2);
    if (num_tags) {
        if (same_prior) {
            AP.compute(0);
            const double* optr=AP.get_offsets();
            output[1]=Rcpp::NumericVector(optr, optr+num_libs);
        } else {
             Rcpp::NumericMatrix current(num_tags, num_libs);
             libptr=&(*current.begin());
             output[1]=current;
        }
    } else {
        if (same_prior) {
             output[1]=Rcpp::NumericVector(num_libs, R_NaReal);
        } else {
             output[1]=Rcpp::NumericMatrix(num_tags, num_libs);
        }
    }

    // Adding the prior values to the existing counts.
    for (int tag=0; tag<num_tags; ++tag) {
        AP.compute(tag);
        const double* pptr=AP.get_priors();

        auto current=outmat.row(tag);
        for (auto& curval : current) { 
            curval += *pptr;
            ++pptr;
        }

        if (!same_prior) {
            const double* optr=AP.get_offsets();
            double* curlibptr=libptr+tag;
            
            for (int lib=0; lib<num_libs; ++lib, curlibptr+=num_tags, ++optr) {
                (*curlibptr)=(*optr);
            }
        }
    }

    output[0] = outmat;
    return output;
    END_RCPP
}

