
 
//Note:  this module was taken from Cam Abrams's code. He took it himself
//from Junichi Tanaka, and added comments.


#ifndef BICUBIC_H
#define BICUBIC_H
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>

#define	X1_NGRIDPOINTS	11
#define X2_NGRIDPOINTS	11
#define X1_NGRIDSQUARES (X1_NGRIDPOINTS-1)
#define X2_NGRIDSQUARES (X2_NGRIDPOINTS-1)

void * xmalloc(size_t n);

/* bcucof:  (Numerical Recipes, 2nd Ed., p 126)  Generates the
 * 16-member coefficient matrix c_(ij), i=1->4, j=1->4.  c_(ij)
 * is used in the "bcuint" function to interpolate a function
 * value at a specified off-lattice point.
 */
void bcucof(double y[], double y1[], double y2[], double y12[],
            double d1, double d2, double c[4][4]);

/* bcuint:  (Numerical Recipes, 2nd Ed., p 127)  A grid square is
 * defined by the four coordinates x1u, x1l, x2u, x2l.  The 
 * 4x4 matrix of coefficients c_(ij) is passed in as c[][].  The
 * point at which an interpolated function value is wished is
 * (x1, x2).  The interpolated function value is returned in *ansy.
 * Its first derivative with respect to x1 is returned in *ansy1.
 * Its first derivative with respect to x2 is returned in *ansy2.
 */
void bcuint(double x1l, double x1u, double x2l, double x2u, 
     double x1, double x2, double *ansy, double *ansy1, double *ansy2,
     double c[4][4]);

/* bicubic_genCoef:  Given the values of the function and its
 * first derivatives and cross derivatives at each grid point, 
 * this function computes the set of all 4x4 coefficient matrices
 * c_(ij) for each grid square. */
void bicubic_genCoef (double y[X1_NGRIDPOINTS][X2_NGRIDPOINTS],
		      double y1[X1_NGRIDPOINTS][X2_NGRIDPOINTS], 
		      double y2[X1_NGRIDPOINTS][X2_NGRIDPOINTS],
		      double y12[X1_NGRIDPOINTS][X2_NGRIDPOINTS], 
		      double (****c)[4][4]);

/* Hcc_genCoef:  Initialize the set of 4x4 coefficient matrices for
 * the Hcc function.  This is the place where the spline knots are
 * specified. */
void Pcc_genCoef (void);
void Pch_genCoef (void);

/* Hcc_bicubicint:  Bicubic interpolation of H_(cc) at point (x1,x2).
 * The value of Hcc at the specified point is returned in *y,
 * d(Hcc)/d(x1) is returned in *y1, d(Hcc)/d(x2) is returned in *y2. */
void Pcc_bicubicint (double x1, double x2,
		     double *y, double *y1, double *y2);

/* Hch_bicubicint:  Bicubic interpolation of H_(ch) at point (x1,x2).
 * The value of Hch at the specified point is returned in *y,
 * d(Hch)/d(x1) is returned in *y1, d(Hch)/d(x2) is returned in *y2. */
void Pch_bicubicint (double x1, double x2,
		     double *y, double *y1, double *y2);

#endif
