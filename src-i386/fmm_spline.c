/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *  Copyright (C) 1998--2001  The R Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/
 */

/* Comment from A. Lun:
 * This is copied straight from 'splines.c' in the 'stats' package of 
 * R's core installation (current as of r72293), with a few modifications:
 *
 * - commented out all functions except for fmm_spline.
 *
 * - commented out the error call in fmm_spline when fewer than 2 points 
 *   are supplied, as that is unnecessary for my code (and fails to compile 
 *   without a definition of EDOM). 
 *
 * - const'ifed the 'x' and 'y' pointers to protect them from modification. 
  */

/*	Spline Interpolation
 *	--------------------
 *	C code to perform spline fitting and interpolation.
 *	There is code here for:
 *	
 *	1. Splines with end-conditions determined by fitting
 *	   cubics in the start and end intervals (Forsythe et al).
 *
 *
 *	Computational Techniques
 *	------------------------
 *	A special LU decomposition for symmetric tridiagonal matrices
 *	is used for all computations, except for periodic splines where
 *	Choleski is more efficient.
 */

/*
 *	Splines a la Forsythe Malcolm and Moler
 *	---------------------------------------
 *	In this case the end-conditions are determined by fitting
 *	cubic polynomials to the first and last 4 points and matching
 *	the third derivitives of the spline at the end-points to the
 *	third derivatives of these cubics at the end-points.
 */

void fmm_spline(int n, const double *x, const double *y, double *b, double *c, double *d)
{
    int nm1, i;
    double t;

    /* Adjustment for 1-based arrays */

    x--; y--; b--; c--; d--;

    if(n < 2) {
//        errno = EDOM;
	return;
    }

    if(n < 3) {
	t = (y[2] - y[1]);
	b[1] = t / (x[2]-x[1]);
	b[2] = b[1];
	c[1] = c[2] = d[1] = d[2] = 0.0;
	return;
    }

    nm1 = n - 1;

    /* Set up tridiagonal system */
    /* b = diagonal, d = offdiagonal, c = right hand side */

    d[1] = x[2] - x[1];
    c[2] = (y[2] - y[1])/d[1];/* = +/- Inf	for x[1]=x[2] -- problem? */
    for(i=2 ; i<n ; i++) {
	d[i] = x[i+1] - x[i];
	b[i] = 2.0 * (d[i-1] + d[i]);
	c[i+1] = (y[i+1] - y[i])/d[i];
	c[i] = c[i+1] - c[i];
    }

    /* End conditions. */
    /* Third derivatives at x[0] and x[n-1] obtained */
    /* from divided differences */

    b[1] = -d[1];
    b[n] = -d[nm1];
    c[1] = c[n] = 0.0;
    if(n > 3) {
	c[1] = c[3]/(x[4]-x[2]) - c[2]/(x[3]-x[1]);
	c[n] = c[nm1]/(x[n] - x[n-2]) - c[n-2]/(x[nm1]-x[n-3]);
	c[1] = c[1]*d[1]*d[1]/(x[4]-x[1]);
	c[n] = -c[n]*d[nm1]*d[nm1]/(x[n]-x[n-3]);
    }

    /* Gaussian elimination */

    for(i=2 ; i<=n ; i++) {
	t = d[i-1]/b[i-1];
	b[i] = b[i] - t*d[i-1];
	c[i] = c[i] - t*c[i-1];
    }

    /* Backward substitution */

    c[n] = c[n]/b[n];
    for(i=nm1 ; i>=1 ; i--)
	c[i] = (c[i]-d[i]*c[i+1])/b[i];

    /* c[i] is now the sigma[i-1] of the text */
    /* Compute polynomial coefficients */

    b[n] = (y[n] - y[n-1])/d[n-1] + d[n-1]*(c[n-1]+ 2.0*c[n]);
    for(i=1 ; i<=nm1 ; i++) {
	b[i] = (y[i+1]-y[i])/d[i] - d[i]*(c[i+1]+2.0*c[i]);
	d[i] = (c[i+1]-c[i])/d[i];
	c[i] = 3.0*c[i];
    }
    c[n] = 3.0*c[n];
    d[n] = d[nm1];
    return;
}
