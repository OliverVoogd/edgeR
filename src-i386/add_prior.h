#ifndef ADD_PRIOR_H
#define ADD_PRIOR_H
#include "objects.h"
#include "utils.h"

class add_prior{
public:
    add_prior(Rcpp::RObject, Rcpp::RObject, bool, bool);
    void compute(int);
    const double* get_priors() const;
    const double* get_offsets() const;

    int get_nrow() const;
    int get_ncol() const; 
    const bool same_across_rows() const;
private:
    compressed_matrix allp, allo;
    const bool logged_in, logged_out;
    int nrow, ncol;

    std::vector<double> adj_prior, adj_libs;
    bool filled;
};

void check_AP_dims(const add_prior&, int, int, const char*);

#endif
