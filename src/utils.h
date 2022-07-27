#ifndef UTILS_H
#define UTILS_H
//#define DEBUG

#ifdef DEBUG
#include <iostream>
#endif

#ifndef USE_FC_LEN_T
#define USE_FC_LEN_T
#endif
#include <Rconfig.h>
#include "R_ext/BLAS.h"
#include "R_ext/Lapack.h"
#ifndef FCONE
#define FCONE
#endif

#include "Rcpp.h"

#include <vector>
#include <cmath>
#include <stdexcept>
#include <sstream>
#include <algorithm>

/* Defining all R-accessible functions. */

extern "C" {

SEXP compute_nbdev(SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP compute_apl (SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP exact_test_by_deviance(SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP fit_levenberg (SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP get_levenberg_start (SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP loess_by_col(SEXP, SEXP, SEXP, SEXP);

SEXP maximize_interpolant(SEXP, SEXP);

SEXP fit_one_group (SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP get_one_way_fitted (SEXP, SEXP, SEXP);

SEXP simple_good_turing (SEXP, SEXP, SEXP);

SEXP add_prior_count (SEXP, SEXP, SEXP);

SEXP calculate_cpm_log (SEXP, SEXP, SEXP);

SEXP calculate_cpm_raw (SEXP, SEXP);

SEXP ave_log_cpm(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP check_poisson_bound (SEXP, SEXP, SEXP);

void processHairpinReads(int *, int *, char**, char**, int*,
		char**, char**, int*, int*, int*, int*, int*, int*,
		int*, int*, int*, int*, int*, int*, int *,
		int *, char**, int*);

}

/* Other utility functions and values */

const double low_value=std::pow(10.0, -10.0), log_low_value=std::log(low_value);

const double LNtwo=std::log(2), one_million=1000000, LNmillion=std::log(one_million);

#endif
