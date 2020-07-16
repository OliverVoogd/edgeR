#include "utils.h"
#include "interpolator.h"

SEXP maximize_interpolant(SEXP spline_pts, SEXP likelihoods) {
    BEGIN_RCPP

    Rcpp::NumericVector spts(spline_pts);
    Rcpp::NumericMatrix ll(likelihoods);
    const int num_pts=spts.size();
    if (num_pts!=ll.ncol()) { 
        throw std::runtime_error("number of columns in likelihood matrix should be equal to number of spline points");
    }
    const int num_tags=ll.nrow();

    interpolator maxinterpol(num_pts);
    std::vector<double> current_ll(num_pts);
    std::vector<double> all_spts(spts.begin(), spts.end()); // making a copy to guarantee contiguousness.

    Rcpp::NumericVector output(num_tags);
    for (int tag=0; tag<num_tags; ++tag) {
        Rcpp::NumericMatrix::Row curll=ll.row(tag);
        std::copy(curll.begin(), curll.end(), current_ll.begin());
        output[tag]=maxinterpol.find_max(all_spts.data(), current_ll.data());
    }

    return(output);    
    END_RCPP
}
