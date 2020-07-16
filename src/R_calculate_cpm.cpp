#include "objects.h"
#include "utils.h"
#include "add_prior.h"

/**** Calculating the CPMs in cpm.default with log=TRUE, but more memory-efficiently. *******/

SEXP calculate_cpm_log (SEXP y, SEXP libsize, SEXP prior) {
    BEGIN_RCPP

    // Checking the inputs.
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

    add_prior AP(prior, libsize, false, true);
    check_AP_dims(AP, num_tags, num_libs, "count");

    // Computing the various statistics.
    for (int tag=0; tag<num_tags; ++tag) {
        AP.compute(tag);
        const double* pptr=AP.get_priors();
        const double* optr=AP.get_offsets();

        auto current_row=outmat.row(tag);
        for (auto& curval : current_row) {     
            curval+=(*pptr);
            if (curval > 0) { 
                curval=std::log(curval)-(*optr)+LNmillion;
                curval/=LNtwo;
            } else {
                curval=R_NaN;
            }

            ++pptr;
            ++optr;
        }
    }
    
    return outmat;
    END_RCPP
}

/**** Calculating the CPMs in cpm.default with log=FALSE, but more memory-efficiently. *******/

SEXP calculate_cpm_raw (SEXP y, SEXP libsize) {
    BEGIN_RCPP

    // Checking the inputs.
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

    compressed_matrix allL(libsize);
    if (allL.get_nrow()!=num_tags || allL.get_ncol()!=num_libs) {
        throw std::runtime_error("dimensions are not consistent between counts and library sizes");
    }

    // Doing the various calculations.
    for (int tag=0; tag<num_tags; ++tag) {
        auto current_row=outmat.row(tag);
        const double* lptr=allL.get_row(tag);

        for (auto& curval : current_row) { 
            if ((*lptr)>0) {
                curval*=one_million/(*lptr);
            } else {
                curval=R_NaN;
            }
            ++lptr;
        }
    }

    return outmat;
    END_RCPP
}



