#include "chemtab.h"

map<short,double> RC_SQ, F_MIN, F_MAX, A, B, MORSE_LAM, MORSE_MU, RE;
map<short,double> BO_ETA, BO_DELTA, MASS, ALPHA, BETA, SPL_RA, SPL_RB;
map<short,double> MO_A, MO_EPS, MO_S, SPL_A, SPL_B, SPL_C;

void ChemDef(){  

  //periodic table of elements
  MASS[6]  = 12.011;    //C
  MASS[18] = 39.948;    //Ar
  MASS[9]  = 18.998403; //F  
  MASS[14] = 28.066;    //Si
  MASS[30] = 4.0026;    //He
  MASS[31] = 20.17;     //Ne
  MASS[17] = 35.453;    //Cl
  
  //carbon-carbon
  RC_SQ[12] = 4.0;
  F_MIN[12] = 1.7;
  F_MAX[12] = 2.0;
  A[12] = 2605.8416;
  B[12] = 1397.073;
  MORSE_LAM[12] = 3.2803;
  MORSE_MU[12] = 2.6888;
  RE[12] = 1.315; //from Tanaka 2000
  SPL_RA[12] = 0.286968;
  SPL_RB[12] = 0.6522;
  SPL_A[12] = -2.862885;
  SPL_B[12] = 7.729378;
  SPL_C[12] = -44.727802;
  MO_A[12] = 0.388874*pow(2*sqrt(6.0),-TWO_THIRDS);
  MO_EPS[12] = 14.39965 * 36;
  MO_S[12] = 544.237667;
  
  //silicon-silicon
  RC_SQ[28] = 9.0;
  F_MIN[28] = 2.7;
  F_MAX[28] = 3.0;
  A[28] = 1830.8;
  B[28] = 471.18;
  MORSE_LAM[28] = 2.4799;
  MORSE_MU[28] = 1.7322;
  RE[28] = 2.35;
  SPL_RA[28] = 0.286968;
  SPL_RB[28] = 0.6522;
  SPL_A[28] = -7.155376;
  SPL_B[28] = 9.502208;
  SPL_C[28] = 237.361562;
  MO_A[28] = 0.388874*pow(2*sqrt(14.0),-TWO_THIRDS);
  MO_EPS[28] = 14.39965 * 14*14;
  MO_S[28] = 296.792932;
  
  //silicon-carbon
  F_MIN[20] = 2.204541;
  F_MAX[20] = 2.509980;
  RC_SQ[20] = pow(F_MAX[20],2);
  A[20] = 1597.311;
  B[20] = 395.1451;
  MORSE_LAM[20] = 2.9839;
  MORSE_MU[20] = 1.972050;
  RE[20] = 1.85;
  SPL_RA[20] = 0.286968;
  SPL_RB[20] = 0.6522;
  SPL_A[20] = -5.9043;
  SPL_B[20] = 8.59831;
  SPL_C[20] = 112.845;
  MO_A[20] = 0.388874*pow(sqrt(14.0)+sqrt(6.0),-TWO_THIRDS);
  MO_EPS[20] = 14.39965 * 14*6;
  MO_S[20] = 292.68;
  
  //silicon-fluorine
  RC_SQ[23] = 4.576262208;
  F_MIN[23] = 1.83922;
  F_MAX[23] = 2.13922;
  SPL_RA[23] = 0.3;
  SPL_RB[23] = 0.8;
  RE[23] = 1.6008;
  MO_A[23] = 0.388874*pow(sqrt(14.0)+sqrt(9.0),-TWO_THIRDS);
  MO_EPS[23] = 14.39965 * 14*9;

// 2003-07-12 -- no Psif function by this date
  // A[23] = 37412.28;
  // B[23] = 925.846;
  // MORSE_MU[23] = 2.7437;
  // MO_S[23] = 1688.552571;
  // SPL_A[23] = -2.1325772404;
  // SPL_B[23] = 8.7909576222;
  // SPL_C[23] = -729.85931698;

//2003-08-09 -- JCP paper submitted Sep 2003
  // A[23] =835.7994893;
  // B[23] =150.2217873;
  // MORSE_MU[23] =1.505136887;
  // SPL_A[23]=-8.873298291;
  // SPL_B[23]=10.34094053;
  // SPL_C[23]=49.67808337;
  // MO_S[23]=-2575.054547;

// 2003-11-07 -- these are the values in JCP paper
  A[23] = 1046.953;
  B[23] = 155.9846;
  MORSE_MU[23] = 1.622332;
  MO_S[23]=-2521.162686;
  SPL_A[23]=-8.654155386;
  SPL_B[23]=10.30020463;
  SPL_C[23]=48.81431439;

  MORSE_LAM[23] = MORSE_MU[23] * 2;

  //silicon-chlorine
  F_MIN[31] = 2.320;
  F_MAX[31] = 2.620;
  RC_SQ[31] = pow(F_MAX[31],2);
  A[31]=	1234.011518;
  B[31]=	142.606135;
  MORSE_MU[31]=	1.411428806;
  MORSE_LAM[31] = 2 * MORSE_MU[31];
  RE[31] = 2.019;
  SPL_RA[31] = 0.3;
  SPL_RB[31] = 0.8;
  MO_A[31] = 0.388874*pow(sqrt(14.0)+sqrt(17.0),-TWO_THIRDS);
  MO_EPS[31] = 14.39965 * 14*17;
  SPL_A[31]=	-9.199255798;
  SPL_B[31]=	11.03777121;
  SPL_C[31]=	89.40890183;
  MO_S[31]=	-4994.987664;

  //chlorine-chlorine
  F_MIN[34] = 2.284;
  F_MAX[34] = 2.584;
  RC_SQ[34] = pow(F_MAX[34],2);
  A[34]=	7248.246586;
  B[34]=	269.7098494;
  MORSE_MU[34]=	2.004384395;
  MORSE_LAM[34] = 2 * MORSE_MU[34];
  RE[34] = 1.9878;
  SPL_RA[34] = 0.3;
  SPL_RB[34] = 0.8;
  MO_A[34] = 0.388874*pow(sqrt(17.0)+sqrt(17.0),-TWO_THIRDS);
  MO_EPS[34] = 14.39965 * 17*17;
  SPL_A[34]=	-7.256044358;
  SPL_B[34]=	10.89298466;
  SPL_C[34]=	131.2995517;
  MO_S[34]=	-4930.327312;

  //fluorine-fluorine FROM TANAKA 2000
  RC_SQ[18] = 4.0;
  F_MIN[18] = 1.7;
  F_MAX[18] = 2.0;
  A[18] = 6379.15;
  B[18] = 2813.8185;
  MORSE_LAM[18] = 4.6407;
  MORSE_MU[18] = 3.9462;
  RE[18] = 1.4119;
  SPL_RA[18] = 0.286968;
  SPL_RB[18] = 0.6522;
  SPL_A[18] = -3.8113458035;
  SPL_B[18] = 8.4167681926;
  SPL_C[18] = -67.291591284;
  MO_A[18] = 0.388874*pow(6.0,-TWO_THIRDS);
  MO_EPS[18] = 14.39965 * 81;
  MO_S[18] = 642.52174315;

  //carbon-fluorine FROM TANAKA 2000
  RC_SQ[15] = 4.0;
  F_MIN[15] = 1.7;
  F_MAX[15] = 2.0;
  A[15] = 1535.1781;
  B[15] = 949.9585;
  MORSE_LAM[15] = 3.1575;
  MORSE_MU[15] = 2.6360;
  RE[15] = 1.2667;
  SPL_RA[15] = 0.286968;
  SPL_RB[15] = 0.6522;
  SPL_A[15] = -5.1521258108;
  SPL_B[15] = 8.1476678001;
  SPL_C[15] = 75.802194262;
  MO_A[15] = 0.388874*pow(3.0+sqrt(6.0),-TWO_THIRDS);
  MO_EPS[15] = 14.39965 * 54.0;
  MO_S[15] = 289.73644378;

  //rest of the Moliere params for inerts.
  MO_A[24] = 0.388874*pow(sqrt(18.0)+sqrt(6.0),-TWO_THIRDS);
  MO_EPS[24] = 14.39965 * 18.0 * 6.0;
  MO_A[27] = 0.388874*pow(sqrt(18.0)+3.0,-TWO_THIRDS);
  MO_EPS[27] = 14.39965 * 18.0 * 9.0;
  MO_A[32] = 0.388874*pow(sqrt(18.0)+sqrt(14.0),-TWO_THIRDS);
  MO_EPS[32] = 14.39965 * 18.0 * 14.0;
  MO_A[35] = 0.388874*pow(sqrt(18.0)+sqrt(17.0),-TWO_THIRDS);
  MO_EPS[35] = 14.39965 * 18.0 * 17.0;
//   MO_A[36] = 0.388874*pow(sqrt(2.0)+sqrt(6.0),-TWO_THIRDS);
//   MO_EPS[36] = 14.39965 * 12.0;
//   MO_A[39] = 0.388874*pow(sqrt(2.0)+3.0,-TWO_THIRDS);
//   MO_EPS[39] = 14.39965 * 18.0;
//   MO_A[44] = 0.388874*pow(sqrt(2.0)+sqrt(14.0),-TWO_THIRDS);
//   MO_EPS[44] = 14.39965 * 28.0;

//   MO_A[37] = 0.388874*pow(sqrt(10.0)+sqrt(6.0),-TWO_THIRDS);
//   MO_EPS[37] = 14.39965 * 60.0;
//   MO_A[40] = 0.388874*pow(sqrt(10.0)+3.0,-TWO_THIRDS);
//   MO_EPS[40] = 14.39965 * 90.0;
//   MO_A[45] = 0.388874*pow(sqrt(10.0)+sqrt(14.0),-TWO_THIRDS);
//   MO_EPS[45] = 14.39965 * 140.0;
}

