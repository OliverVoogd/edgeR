#include "utils.h"
#include "glm.h"
#include "objects.h"

SEXP exact_test_by_deviance(SEXP sums_1, SEXP sums_2, SEXP n_1, SEXP n_2, SEXP disp) {
    BEGIN_RCPP
    Rcpp::IntegerVector S1(sums_1), S2(sums_2), dispersions(disp);
    const int num_tags=S1.size();
    if (num_tags!=S2.size() || num_tags!=dispersions.size()) { 
        throw std::runtime_error("lengths of input vectors do not match");
    }

    const int n1=check_integer_scalar(n_1, "number of libraries");
    const int n2=check_integer_scalar(n_2, "number of libraries");
    if (n1<=0 || n2 <=0) { 
        throw std::runtime_error("number of libraries must be positive for each condition");
    }
	const int nlibs = n1+n2;

    Rcpp::NumericVector output(num_tags);
    for (int i=0; i<num_tags; ++i) {
        const int& s1=S1[i];
 	    const int& s2=S2[i];
        const int stotal=s1+s2;

		// Computing current means and sizes for each library (probability is the same).
		const double mu = stotal/nlibs;
		const double mu1=mu*n1, mu2=mu*n2, r1=n1/dispersions[i], r2=n2/dispersions[i];
        const double p = r1/(r1+mu1);

		/* The aim is to sum conditional probabilities for all partitions of the total sum with deviances 
 		 * greater than that observed for the current partition. We start computing from the extremes
 		 * in both cases.
 		 */
		const double phi1=1/r1, phi2=1/r2;
		const double obsdev=compute_unit_nb_deviance(s1, mu1, phi1)+compute_unit_nb_deviance(s2, mu2, phi2);
		double& currentp=output[i];
		
		// Going from the left.	
		int j=0;
		while (j <= stotal) {
			if (obsdev <= compute_unit_nb_deviance(j, mu1, phi1)+compute_unit_nb_deviance(stotal-j, mu2, phi2)) { 
				currentp+=R::dnbinom(j, r1, p, 0) * R::dnbinom(stotal-j, r2, p, 0);
			} else { break; }
			++j;
		}

		// Going from the right, or what's left of it.
		for (int k=0; k<=stotal-j; ++k) {
			if (obsdev <= compute_unit_nb_deviance(k, mu2, phi2)+compute_unit_nb_deviance(stotal-k, mu1, phi1)) { 
				currentp+=R::dnbinom(k, r2, p, 0) * R::dnbinom(stotal-k, r1, p, 0);
			} else { break; }
		}

    	const double totalr=r1+r2;
	    currentp /= R::dnbinom(stotal, totalr, totalr/(totalr+mu1+mu2), 0);
    }
    
    return output;
    END_RCPP
}
