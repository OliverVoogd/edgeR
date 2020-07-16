#include "glm.h"
#include "objects.h"

/*** Function to compute the fitted values without a lot of temporary matrices. ***/

SEXP get_one_way_fitted (SEXP beta, SEXP offset, SEXP groups) { 
    BEGIN_RCPP
    Rcpp::NumericMatrix Beta(beta);
    const int num_tags=Beta.nrow();
    const int num_groups=Beta.ncol();
  
    Rcpp::IntegerVector Groups(groups);
    const int num_libs=Groups.size();
    if (*std::min_element(Groups.begin(), Groups.end()) < 0) {
        throw std::runtime_error("smallest value of group vector should be non-negative");
    }
    if (*std::max_element(Groups.begin(), Groups.end()) >= num_groups) {
        throw std::runtime_error("largest value of group vector should be less than the number of groups");
    }
 
    compressed_matrix allo=check_CM_dims(offset, num_tags, num_libs, "offset", "count");

    Rcpp::NumericMatrix output(num_tags, num_libs);
    std::vector<double> betabyrow(num_libs);
    for (int tag=0; tag<num_tags; ++tag) {
        auto curbeta=Beta.row(tag);
        std::copy(curbeta.begin(), curbeta.end(), betabyrow.begin());
        const double* optr=allo.get_row(tag);
       
        auto gIt=Groups.begin(); 
        auto curout=output.row(tag);
        auto coIt=curout.begin();
        for (int lib=0; lib<num_libs; ++lib, ++coIt, ++gIt) {
            (*coIt)=std::exp(optr[lib] + betabyrow[*gIt]);
        }
    }

    return output;
    END_RCPP
}

