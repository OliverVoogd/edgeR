#include "R_ext/Rdynload.h"
#include "R_ext/Visibility.h"
#include "utils.h"

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

extern "C" {

static const R_CallMethodDef all_call_entries[] = {
	CALLDEF(compute_nbdev, 5),
	CALLDEF(compute_apl, 6),
	CALLDEF(exact_test_by_deviance, 5), 
	CALLDEF(loess_by_col, 4),
	CALLDEF(maximize_interpolant, 2),

    CALLDEF(fit_levenberg, 8),
	CALLDEF(get_levenberg_start, 6),
	CALLDEF(fit_one_group, 7),
	CALLDEF(get_one_way_fitted, 3),
	CALLDEF(simple_good_turing, 3),

    CALLDEF(add_prior_count, 3),
    CALLDEF(calculate_cpm_log, 3),
    CALLDEF(calculate_cpm_raw, 2),
    CALLDEF(ave_log_cpm, 7),

    CALLDEF(check_poisson_bound, 3),
	{NULL, NULL, 0}
};

R_CMethodDef all_c_entries[] = {
    {"processHairpinReads", (DL_FUNC) &processHairpinReads, 20},
    {NULL, NULL, 0}
  };

void attribute_visible R_init_edgeR(DllInfo *dll) {
	R_registerRoutines(dll, all_c_entries, all_call_entries, NULL, NULL);
	R_useDynamicSymbols(dll, FALSE);
	R_forceSymbols(dll, TRUE);
}

}
