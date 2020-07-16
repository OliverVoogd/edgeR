#include "interpolator.h"

struct solution {
	double sol1, sol2;
	bool solvable;
};

solution quad_solver (const double& a, const double& b, const double& c) {
	double back=b*b-4*a*c;
	solution cur_sol;
	if (back<0) {
		cur_sol.solvable=false;
		return(cur_sol);
	}
	double front=-b/(2*a);
	back=std::sqrt(back)/(2*a);
	cur_sol.sol1=front-back;
	cur_sol.sol2=front+back;
	cur_sol.solvable=true;
	return(cur_sol);
}

/************************************
 *  
 * It fits the spline and grabs the coefficients of each segment.
 * It then calculates the maxima at the segments neighbouring
 * the maximally highest point. This avoids NR optimisation 
 * as well as the need to call R's splinefun's from within C.
 * 
 ***********************************/

interpolator::interpolator(const int& n) : npts(n), b(npts), c(npts), d(npts) {
	if (npts<2) { throw std::runtime_error("must have at least two points for interpolation"); }
    return;
}

double interpolator::find_max (const double*x, const double* y) {
    double maxed=-1;
	int maxed_at=-1;
	for (int i=0; i<npts; ++i) {
	// Getting a good initial guess for the MLE.
	    if (maxed_at < 0 || y[i] > maxed) {
           	maxed=y[i];
           	maxed_at=i;
 	   	}
	}
    double x_max=x[maxed_at];
    fmm_spline(npts, x, y, b.data(), c.data(), d.data());

	// First we have a look at the segment on the left and see if it contains the maximum.
    if (maxed_at>0) {
        const double& ld=d[maxed_at-1];
        const double& lc=c[maxed_at-1];
        const double& lb=b[maxed_at-1];
        solution sol_left=quad_solver(3*ld, 2*lc, lb);
        if (sol_left.solvable) {
            /* Using the solution with the maximum (not minimum). If the curve is mostly increasing, the 
             * maximal point is located at the smaller solution (i.e. sol1 for a>0). If the curve is mostly
             * decreasing, the maximal is located at the larger solution (i.e., sol1 for a<0).
             */
            const double& chosen_sol=sol_left.sol1;

            /* The spline coefficients are designed such that 'x' in 'y + b*x + c*x^2 + d*x^3' is
             * equal to 'x_t - x_l' where 'x_l' is the left limit of that spline segment and 'x_t'
             * is where you want to get an interpolated value. This is necessary in 'splinefun' to 
             * ensure that you get 'y' (i.e. the original data point) when 'x=0'. For our purposes, 
             * the actual MLE corresponds to 'x_t' and is equal to 'solution + x_0'.
             */
            if (chosen_sol > 0 && chosen_sol < x[maxed_at]-x[maxed_at-1]) {
                double temp=((ld*chosen_sol+lc)*chosen_sol+lb)*chosen_sol+y[maxed_at-1];
                if (temp > maxed) {
                    maxed=temp;
                    x_max=chosen_sol+x[maxed_at-1];
                }
            }
        }
    }
    
	// Repeating for the segment on the right.
    if (maxed_at<npts-1) {
        const double& rd=d[maxed_at];
        const double& rc=c[maxed_at];
        const double& rb=b[maxed_at];
        solution sol_right=quad_solver(3*rd, 2*rc, rb);
        if (sol_right.solvable) {
            const double& chosen_sol=sol_right.sol1; // see arguments above.
            if (chosen_sol > 0 && chosen_sol < x[maxed_at+1]-x[maxed_at]) {
                double temp=((rd*chosen_sol+rc)*chosen_sol+rb)*chosen_sol+y[maxed_at];
                if (temp>maxed) {
                    maxed=temp;
                    x_max=chosen_sol+x[maxed_at];
                }
            }
        }
    }

	return x_max;
}


