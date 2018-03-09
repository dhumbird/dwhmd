//  Note:  Junichi Tanaka did all of the hard work for this module.
//  I stole it from Cam, who borrowed it from him.

#include "tricubic.h"
#include <iostream>
static const int mp = (X1_NGRIDPOINTS);
static const int np = (X2_NGRIDPOINTS);
static const int ms = (X1_NGRIDSQUARES);
static const int ns = (X2_NGRIDSQUARES);
static const int lp = (X3_NGRIDPOINTS);
static const int ls = (X3_NGRIDSQUARES);

int TRICUBIC_DIAG_=0;
short fcc_opt_=0;

static  double (****cTCC)[4][4][4];
static  double (****cCC)[4][4][4];
static  double (****cHH)[4][4][4];
static  double (****cCH)[4][4][4];
/* A ptr to...__||||    |_______|
 *	         |||         |_________________________ 
 *               |||                                   |
 *		 \|/				       |
 *		  |____an "m"x"n"x"l" 3D array of...   |___ 4x4x4 3D arrays.
 *
 *  Here, we state that the 3D grid of points ("knots") which define the
 *  function we are interpolating is "m+1"x"n+1"x"l+1" and has, therefore, 
 *  m x n x l 8-membered grid cubes.  Each one of these grid cubes 
 *  has associated with it a 4x4x4 matrix of coefficients, c_(ijk), 
 *  which are functions of the function values and its derivatives
 *  at the eight gridpoints in that cube.  This coefficient matrix
 *  c_(ijk) is used in the tricubic interpolation routine.  We wish to
 *  initially compute *all* 4x4x4 matrices (c_(ijk))_(mnl) initially, and
 *  store them for use by the interpolation routine.  In order to store
 *  the m x n x l 4x4x4 arrays, we declare "c" as a "(pointer to a) 
 *  3D-matrix of 3D-matrices",  or (*)(***c)[4][4][4], and dynamically
 *  allocate the "m" superrows of "n" columns of "l" rows each using the
 *  xmalloc function adapted from the bicubic analog in Numerical Recipes.
 *
 */

void tcucof(double y[8], double y1[8], double y2[8], double y3[8],
	    double y12[8], double y23[8], double y31[8], double y123[8],
            double d1, double d2, double d3, double c[4][4][4]){
  /* Declare the 4,096 tabulated weighting factors in a [64]x[64]
   * array. */
#include "tricof.h"
  
  int l, k, j, i;
  double xx, cl[64], x[64], d1d2, d2d3, d3d1, d1d2d3;
  
  d1d2 = d1*d2;
  d2d3 = d2*d3;
  d3d1 = d3*d1;
  d1d2d3 = d1d2*d3;
  
  /* Build the temporary vector */
  for (i = 0; i < 8; i++){
    x[i]	=   y[i];
    x[i+8]	=   y1[i]*d1;
    x[i+16]	=   y2[i]*d2;
    x[i+24]	=   y3[i]*d3;
    x[i+32]	=   y12[i]*d1d2;
    x[i+40]	=   y23[i]*d2d3;
    x[i+48]	=   y31[i]*d3d1;
    x[i+56]	=   y123[i]*d1d2d3;
  }
  
  /* Matrix multiply the weighting factors with the temporary vector,
   * and store the result in a second vector. */
  for (i = 0; i < 64; i++){
    xx = 0.0;
    for (k = 0; k < 64; k++) xx += twt[i][k]*x[k];
    cl[i] = xx;
  }
    
  /* Unpack the second temporary vector into the 4x4x4 coefficient matrix. */
  l = 0;
  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++)
      for (k = 0; k < 4; k++)
	c[i][j][k] = cl[l++];
}

void tcuint(double x1l, double x1u, double x2l, double x2u, 
	    double x3l, double x3u, 
	    double x1, double x2, double x3, 
	    double *ansy, double *ansy1, double *ansy2, double *ansy3, 
	    double c[4][4][4]){
  int i, j, k;
  double t, u, v, d1, d2, d3;
  double q1, q2, q3;
  d1 = x1u - x1l;
  d2 = x2u - x2l;
  d3 = x3u - x3l;

  if (!d1 || !d2 || !d3){
    fprintf(stderr, "3D Cubic Interpolation: Bad gridpoint coords:\n");
    fprintf(stderr, "\tx1u(%.5e) x1l(%.5e) x2u(%.5e) x2l(%.5e) x3u(%.5e) x3l(%.5e)\n", 
	    x1u, x1l, x2u, x2l, x3u, x3l);
    fprintf(stderr, "Program exits\n");
    exit(0);
  }
  t = (x1 - x1l)/d1;
  u = (x2 - x2l)/d2;
  v = (x3 - x3l)/d3;
  *ansy = (*ansy1) = (*ansy2) = (*ansy3) = 0.0;
  
  for (i = 3; i >= 0; i--){
    q1=(((c[i][3][3]*v + c[i][3][2])*v + c[i][3][1])*v + c[i][3][0])*u*u;
    q2=(((c[i][2][3]*v + c[i][2][2])*v + c[i][2][1])*v + c[i][2][0])*u;
    q3=(((c[i][1][3]*v + c[i][1][2])*v + c[i][1][1])*v + c[i][1][0]);
    *ansy = t*(*ansy)
      + (q1 + q2 + q3)*u
      + (((c[i][0][3]*v + c[i][0][2])*v + c[i][0][1])*v + c[i][0][0]);
    *ansy1 = u*(*ansy1)
      + ((((c[3][i][3]*v + c[3][i][2])*v + c[3][i][1])*v + c[3][i][0])*3*t
	 + (((c[2][i][3]*v + c[2][i][2])*v + c[2][i][1])*v + c[2][i][0])*2)*t
      + (((c[1][i][3]*v + c[1][i][2])*v + c[1][i][1])*v + c[1][i][0]);
    *ansy2 = t*(*ansy2) + 3*q1 + 2*q2 + q3;
    *ansy3 = t*(*ansy3)
      + ((((c[i][3][3]*u + c[i][2][3])*u + c[i][1][3])*u + c[i][0][3])*3*v
	 + (((c[i][3][2]*u + c[i][2][2])*u + c[i][1][2])*u + c[i][0][2])*2)*v
      + (((c[i][3][1]*u + c[i][2][1])*u + c[i][1][1])*u + c[i][0][1]);
  }
  *ansy1 /= d1;
  *ansy2 /= d2;
  *ansy3 /= d3;    
}

void tricubic_genCoef 
(double y[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS],
 double y1[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS], 
 double y2[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS],
 double y3[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS],
 double y12[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS], 
 double y23[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS], 
 double y31[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS], 
 double y123[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS], 
 double (*****c)[4][4][4]){
  int i, j, k, s, ii, jj, kk;
  double z[8], z1[8], z2[8], z3[8], z12[8], z23[8], z31[8], z123[8];
  int ip[][3] = {{0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,1},{1,0,1},
		 {1,1,1},{0,1,1}};
  
    /* The array ip[][3] determines the order of visitation of the 8 vertices
     * of the cube:
     * 
     * 
     *              (8)------(7)
     *              /|       /|
     *             / |      / |		 Directions:
     *		(4)------(3)  |
     *		 |  (5)---|--(6)               /\    __ 
     *           |  /     |  /                 ||    //|
     *           | /      | /                  ||   //
     *          (1)------(2)             1===>  2  3 
     * 
     */

  /* Initialize *c as a matrix of matrices */
  (*c) = (double (****)[4][4][4])xmalloc(ls*sizeof(double (***)[4][4][4]));
  for (i = 0; i < lp; i++){
    (*c)[i] = (double (***)[4][4][4])xmalloc(ms*sizeof(double (**)[4][4][4]));
    for (j = 0; j < mp; j++){
      (*c)[i][j] = (double (**)[4][4][4])xmalloc(ns*sizeof(double (*)[4][4][4]));
      for (k = 0; k < np; k++)
	(*c)[i][j][k] = (double (*)[4][4][4])xmalloc(64*sizeof(double));
    }
  }
  
    for (i = 0; i < ls; i++){
      for (j = 0; j < ms; j++){
	for (k = 0; k < ns; k++){
	  for (s = 0; s < 8; s++){
	    ii = i + ip[s][0];
	    jj = j + ip[s][1];
	    kk = k + ip[s][2];
	    if (ii>=X1_NGRIDPOINTS) {printf("error big_i\n");exit(0);}
	    if (jj>=X2_NGRIDPOINTS) {printf("error big_j\n");exit(0);}
	    if (kk>=X3_NGRIDPOINTS) {printf("error big_k\n");exit(0);}
	    if (ii<0) {printf("error little_i\n");exit(0);}
	    if (jj<0) {printf("error little_j\n");exit(0);}
	    if (kk<0) {printf("error little_k\n");exit(0);}
	    z   [s] = y   [ii][jj][kk];
	    z1  [s] = y1  [ii][jj][kk];
	    z2  [s] = y2  [ii][jj][kk];
	    z3  [s] = y3  [ii][jj][kk];
	    z12 [s] = y12 [ii][jj][kk];
	    z23 [s] = y23 [ii][jj][kk];
	    z31 [s] = y31 [ii][jj][kk];
	    z123[s] = y123[ii][jj][kk];
	  }
	  /* Use tcucof() to compute the 4x4x4 coefficient matrix c_(ijk) */
	  tcucof(z,z1,z2,z3,z12,z23,z31,z123,1.0,1.0,1.0,*((*c)[i][j][k]));
	}
      }
    }
}

void Fcc_genCoef (void){
  int i, j, k;
  double fcc[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS];
  double fcc1[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS];
  double fcc2[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS];
  double fcc3[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS];
  double fcc12[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS];
  double fcc23[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS];
  double fcc31[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS];
  double fcc123[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS];
  
  for (i = 0; i < lp; i++)
    for (j = 0; j < mp; j++)
      for (k = 0; k < np; k++)
	fcc[i][j][k] = 
	  fcc1[i][j][k] = fcc2[i][j][k] = fcc3[i][j][k] =
	  fcc12[i][j][k] = fcc23[i][j][k] = fcc31[i][j][k] =
	  fcc123[i][j][k] = 0.0;

    /*
     * fcc[a][b][c]:  a = "coordination" of c_i excluding j, 
     *		      b = "coordination" of c_j excluding j, 
     *		      c = Nij_conj.
     */


  fcc[0][0][3] = 0.0099172158;
  fcc[0][0][4] = 0.0099172158;
  fcc[0][0][5] = 0.0099172158;
  fcc[0][0][6] = 0.0099172158;
  fcc[0][0][7] = 0.0099172158;
  fcc[0][0][8] = 0.0099172158;
  fcc[0][0][9] = 0.0099172158;
  fcc[0][1][1] = 0.04338699;
  fcc[0][1][2] = 0.0099172158;
  fcc[0][1][3] = 0.0099172158;
  fcc[0][1][4] = 0.0099172158;
  fcc[0][1][5] = 0.0099172158;
  fcc[0][1][6] = 0.0099172158;
  fcc[0][1][7] = 0.0099172158;
  fcc[0][1][8] = 0.0099172158;
  fcc[0][1][9] = 0.0099172158;
  fcc[0][2][1] = 0.0493976637;
  fcc[0][2][2] =-0.011942669;
  fcc[0][2][3] = 0.0099172158;
  fcc[0][2][4] = 0.0099172158;
  fcc[0][2][5] = 0.0099172158;
  fcc[0][2][6] = 0.0099172158;
  fcc[0][2][7] = 0.0099172158;
  fcc[0][2][8] = 0.0099172158;
  fcc[0][2][9] = 0.0099172158;
  fcc[0][3][1] =-0.119798935;
  fcc[0][3][2] =-0.119798935;
  fcc[0][3][3] = 0.0099172158;
  fcc[0][3][4] = 0.0099172158;
  fcc[0][3][5] = 0.0099172158;
  fcc[0][3][6] = 0.0099172158;
  fcc[0][3][7] = 0.0099172158;
  fcc[0][3][8] = 0.0099172158;
  fcc[0][3][9] = 0.0099172158;
  fcc[1][1][1] = 0.105;
  fcc[1][1][2] =-0.0041775;
  fcc[1][1][3] =-0.0160856;
  fcc[1][1][4] =-0.0160856;
  fcc[1][1][5] =-0.0160856;
  fcc[1][1][6] =-0.0160856;
  fcc[1][1][7] =-0.0160856;
  fcc[1][1][8] =-0.0160856;
  fcc[1][1][9] =-0.0160856;
  fcc[1][2][1] = 0.0096495698;
  fcc[1][2][2] = 0.03;
  fcc[1][2][3] =-0.02;
  fcc[1][2][4] =-0.0233778774;
  fcc[1][2][5] =-0.0267557548;
  fcc[1][2][6] =-0.030133632;
  fcc[1][2][7] =-0.030133632;
  fcc[1][2][8] =-0.030133632;
  fcc[1][2][9] =-0.030133632;
  fcc[1][3][2] =-0.124836752;
  fcc[1][3][3] =-0.124836752;
  fcc[1][3][4] =-0.124836752;
  fcc[1][3][5] =-0.124836752;
  fcc[1][3][6] =-0.124836752;
  fcc[1][3][7] =-0.124836752;
  fcc[1][3][8] =-0.124836752;
  fcc[1][3][9] =-0.124836752;
  fcc[2][2][1] = 0.09444957;
  fcc[2][2][2] = 0.022;
  fcc[2][2][3] = 0.03970587;
  fcc[2][2][4] = 0.03308822;
  fcc[2][2][5] = 0.02647058;
  fcc[2][2][6] = 0.01985293;
  fcc[2][2][7] = 0.01323529;
  fcc[2][2][8] = 0.00661764;
  fcc[2][3][1] =-0.044709383;
  fcc[2][3][2] =-0.044709383;
  fcc[2][3][3] =-0.044709383;
  fcc[2][3][4] =-0.044709383;
  fcc[2][3][5] =-0.044709383;
  fcc[2][3][6] =-0.044709383;
  fcc[2][3][7] =-0.044709383;
  fcc[2][3][8] =-0.044709383;
  fcc[2][3][9] =-0.044709383;

  
  fcc1[1][3][2] = 0.0375447763847205;
  fcc1[2][1][1] =-0.0525002674187363;
  fcc1[2][1][5] =-0.05437559;
  fcc1[2][1][6] =-0.05437559;
  fcc1[2][1][7] =-0.05437559;
  fcc1[2][1][8] =-0.05437559;
  fcc1[2][1][9] =-0.05437559;
  fcc1[2][3][2] = 0.0624184;
  fcc1[2][3][3] = 0.0624184;
  fcc1[2][3][4] = 0.0624184;
  fcc1[2][3][5] = 0.0624184;
  fcc1[2][3][6] = 0.0624184;
  fcc1[2][3][7] = 0.0624184;
  fcc1[2][3][8] = 0.0624184;
  fcc1[2][3][9] = 0.0624184;

  fcc3[1][1][2]=-0.0605430516924343;
  fcc3[1][2][4]=-0.020044544;
  fcc3[1][2][5]=-0.020044544;
  fcc3[2][2][4]=-0.0066176451;
  fcc3[2][2][5]=-0.0066176451;
  fcc3[2][2][6]=-0.0066176451;
  fcc3[2][2][7]=-0.0066176451;
  fcc3[2][2][8]=-0.0066176451;

  /* Enforce the symmetry of Fcc */
  for (i=0;i<lp;i++)
    for (j=0;j<mp;j++)
      for(k=0;k<np;k++){
	if (i<j){
 	  fcc[j][i][k]=fcc[i][j][k];
	  fcc3[j][i][k]=fcc3[i][j][k];
	}
	fcc2[j][i][k]=fcc1[i][j][k];
      }
  tricubic_genCoef(fcc, fcc1, fcc2, fcc3, fcc12, fcc23, fcc31, fcc123, &cCC);
//   double y, y1,y2,y3;
//   Fcc_tricubicint(0.3, 1, 2.1, &y, &y1, &y2, &y3);
//   cerr<<y<<" "<<y1<<" "<<y2<<" "<<y3<<endl;
}

void Fcc_tricubicint (double x1, double x2, double x3,
		      double *y, double *y1, double *y2, double *y3){
  int i, j, k;
  double x1max = ls-1e-10, x2max = ms-1e-10, x3max = ns-1e-10;
  double t=*y;
  x1 = (x1 > x1max ? x1max : x1);
  x2 = (x2 > x2max ? x2max : x2);
  x3 = (x3 > x3max ? x3max : x3);
  i = (int)x1;
  j = (int)x2;
  k = (int)x3;
  int ii=(int)*y1;
  int jj=(int)*y2;
  int kk=(int)*y3;
  if (i>=3){i=3; x1=3;}
  if (j>=3){j=3; x2=3;}
  if (k>=9){k=9; x3=9;}
  /* Perform the interpolation on a unit cube:  This makes
   * "t", "u", and "v" 1.0 in the tcuint routine and thus saves time
   * by simplifying multiplication operations. */
  tcuint(0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 
	 x1-i, x2-j, x3-k, 
	 y, y1, y2, y3, *(cCC[i][j][k]));
//    if (*y1+*y2+*y3!=0) 
//      printf("%11.7f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f\n", t, x1, x2, x3, *y, *y1, *y2, *y3);
}
//**************************************************************************
void Fhh_genCoef (void){
  int i, j, k;
  double fhh[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS];
  double fhh1[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS];
  double fhh2[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS];
  double fhh3[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS];
  double fhh12[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS];
  double fhh23[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS];
  double fhh31[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS];
  double fhh123[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS];
  
  for (i = 0; i < lp; i++)
    for (j = 0; j < mp; j++)
      for (k = 0; k < np; k++)
	fhh[i][j][k] = 
	  fhh1[i][j][k] = fhh2[i][j][k] = fhh3[i][j][k] =
	  fhh12[i][j][k] = fhh23[i][j][k] = fhh31[i][j][k] =
	  fhh123[i][j][k] = 0.0;

  fhh[1][1][1] = 0.249831916;
  tricubic_genCoef(fhh, fhh1, fhh2, fhh3, fhh12, fhh23, fhh31, fhh123, &cHH);
//   double y, y1,y2,y3;
//   Fhh_tricubicint(1, 1, 1.1, &y, &y1, &y2, &y3);
//   cerr<<y<<" "<<y1<<" "<<y2<<" "<<y3<<endl;
}

void Fhh_tricubicint (double x1, double x2, double x3,
		      double *y, double *y1, double *y2, double *y3){
  int i, j, k;
  double x1max = ls-1e-10, x2max = ms-1e-10, x3max = ns-1e-10;
  
  x1 = (x1 > x1max ? x1max : x1);
  x2 = (x2 > x2max ? x2max : x2);
  x3 = (x3 > x3max ? x3max : x3);
  i = (int)x1;
  j = (int)x2;
  k = (int)x3;

  /* Perform the interpolation on a unit cube:  This makes
   * "t", "u", and "v" 1.0 in the tcuint routine and thus saves time
   * by simplifying multiplication operations. */
  tcuint(0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 
	 x1-i, x2-j, x3-k, 
	 y, y1, y2, y3, *(cHH[i][j][k]));
}
//*********************************************************************
void Fch_genCoef (void){
  int i, j, k;
  double fch[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS];
  double fch1[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS];
  double fch2[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS];
  double fch3[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS];
  double fch12[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS];
  double fch23[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS];
  double fch31[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS];
  double fch123[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS];
  
  for (i = 0; i < lp; i++)
    for (j = 0; j < mp; j++)
      for (k = 0; k < np; k++)
	fch[i][j][k] = 
	  fch1[i][j][k] = fch2[i][j][k] = fch3[i][j][k] =
	  fch12[i][j][k] = fch23[i][j][k] = fch31[i][j][k] =
	  fch123[i][j][k] = 0.0;

  fch[0][2][5]=-0.00904779;
  fch[0][2][6]=-0.00904779;
  fch[0][2][7]=-0.00904779;
  fch[0][2][8]=-0.00904779;
  fch[0][2][9]=-0.00904779;
  fch[1][1][1]=-0.10;
  fch[1][1][2]=-0.1;
  fch[1][1][3]=-0.6;
  fch[1][1][4]=-0.1;
  fch[1][1][5]=-0.1;
  fch[1][1][6]=-0.1;
  fch[1][1][7]=-0.1;
  fch[1][1][8]=-0.1;
  fch[1][1][9]=-0.1;
  fch[1][2][2]=-0.5;
  fch[1][2][3]=-0.5;
  fch[1][3][1]=-0.2;
  fch[1][3][2]=-0.25;
  fch[1][3][3]=-0.25;
  fch[1][3][4]=-0.2;
  fch[1][3][5]=-0.2;
  fch[1][3][6]=-0.2;
  fch[1][3][7]=-0.2;
  fch[1][3][8]=-0.2;
  fch[1][3][9]=-0.2;

  /* Enforce the symmetry of Fch */
  for (i=0;i<lp;i++)
    for (j=0;j<mp;j++)
      for(k=0;k<np;k++)
	if (i<j) fch[j][i][k]=fch[i][j][k];
  tricubic_genCoef(fch, fch1, fch2, fch3, fch12, fch23, fch31, fch123, &cCH); 
}

void Fch_tricubicint (double x1, double x2, double x3,
		      double *y, double *y1, double *y2, double *y3){
  int i, j, k;
  double x1max = ls-1e-10, x2max = ms-1e-10, x3max = ns-1e-10;
  
  x1 = (x1 > x1max ? x1max : x1);
  x2 = (x2 > x2max ? x2max : x2);
  x3 = (x3 > x3max ? x3max : x3);
  i = (int)x1;
  j = (int)x2;
  k = (int)x3;
  if (i>=3){i=3; x1=3;}
  if (j>=3){j=3; x2=3;}
  if (k>=9){k=9; x3=9;}

  /* Perform the interpolation on a unit cube:  This makes
   * "t", "u", and "v" 1.0 in the tcuint routine and thus saves time
   * by simplifying multiplication operations. */
  tcuint(0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 
	 x1-i, x2-j, x3-k, 
	 y, y1, y2, y3, *(cCH[i][j][k]));
}
//*********************************************************************
void Tcc_genCoef (void){
  int i, j, k;
  double tcc[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS];
  double tcc1[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS];
  double tcc2[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS];
  double tcc3[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS];
  double tcc12[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS];
  double tcc23[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS];
  double tcc31[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS];
  double tcc123[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS];
  
  for (i = 0; i < lp; i++)
    for (j = 0; j < mp; j++)
      for (k = 0; k < np; k++)
	tcc[i][j][k] = 
	  tcc1[i][j][k] = tcc2[i][j][k] = tcc3[i][j][k] =
	  tcc12[i][j][k] = tcc23[i][j][k] = tcc31[i][j][k] =
	  tcc123[i][j][k] = 0.0;

  /* tcc[a][b][c]:  a = "coordination" of c_i excluding j, 
   *		      b = "coordination" of c_j excluding j, 
   *		      c = Nij_conj.
   */
  tcc1[1][2][2]=tcc2[2][1][2]=-0.004048375;
  tcc[2][2][1]=-0.0702801;
  for(k=2;k<np;k++) tcc[2][2][k]=-0.00809675;
  
  tricubic_genCoef(tcc, tcc1, tcc2, tcc3, tcc12, tcc23, tcc31, tcc123, &cTCC); 
}

void Tcc_tricubicint (double x1, double x2, double x3,
		      double *y, double *y1, double *y2, double *y3){
  int i, j, k;
  double x1max = ls-1e-10, x2max = ms-1e-10, x3max = ns-1e-10;
  
  x1 = (x1 > x1max ? x1max : x1);
  x2 = (x2 > x2max ? x2max : x2);
  x3 = (x3 > x3max ? x3max : x3);
  i = (int)x1;
  j = (int)x2;
  k = (int)x3;
  if (i>=3){i=3; x1=3;}
  if (j>=3){j=3; x2=3;}
  if (k>=9){k=9; x3=9;}

  /* Perform the interpolation on a unit cube:  This makes
   * "t", "u", and "v" 1.0 in the tcuint routine and thus saves time
   * by simplifying multiplication operations. */
  tcuint(0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 
	 x1-i, x2-j, x3-k, 
	 y, y1, y2, y3, *(cTCC[i][j][k]));
}
