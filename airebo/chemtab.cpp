#include "chemtab.h"

map<short,double> MASS, RC_SQ, REBO_F_MIN, REBO_F_MAX, Q, A, 
ALPHA, B1, B2, B3, BETA1, REBO_W_MIN, REBO_W_MAX,
LJ_F_MIN, LJ_F_MAX, B_MIN, B_MAX, LJ_SIGMA, LJ_EPS,
BETA2, BETA3, LJ_RCSQ, V_LJ_RC;

void ChemDef(){  

  //type 0: carbon-carbon
  RC_SQ[0] = 4.0;
  REBO_F_MIN[0] = REBO_W_MIN[0] = 1.7;
  REBO_F_MAX[0] = REBO_W_MAX[0] = 2.0;
  Q[0] = 0.3134602960833;
  A[0] = 10953.544162170;
  ALPHA[0] = 4.7465390606595;
  B1[0] = 12388.79197798;
  B2[0] = 17.56740646509;
  B3[0] = 30.71493208065;
  BETA1[0] = 4.7204523127;
  BETA2[0] = 1.4332132499;
  BETA3[0] = 1.3826912506;
  LJ_SIGMA[0] = 3.4;
  LJ_EPS[0] = 4* 0.00284;
  LJ_F_MIN[0] = LJ_SIGMA[0];
  LJ_F_MAX[0] = TWO_POW_ONE_SIX*LJ_SIGMA[0];
  B_MIN[0] = 0.77;
  B_MAX[0] = 0.81;
  LJ_RCSQ[0] = pow(3 * LJ_SIGMA[0],2);
  V_LJ_RC[0] = LJ_EPS[0] * (pow(ONE_THIRD, 12)-pow(ONE_THIRD,6));
 
  //type 1: hydrogen-hydrogen
  RC_SQ[1] = 2.89;
  REBO_F_MIN[1] = REBO_W_MIN[1] = 1.1;
  REBO_F_MAX[1] = REBO_W_MAX[1] = 1.7;
  Q[1] = 0.370471487045;
  A[1] = 32.817355747;
  ALPHA[1] = 3.536298648;
  B1[1] = 29.632593;
  B2[1] = 0.0;
  B3[1] = 0.0;
  BETA1[1] = 1.71589217;
  BETA2[1] = 0.0;
  BETA3[1] = 0.0;
  LJ_SIGMA[1] = 2.65;
  LJ_EPS[1] = 4 * 0.0015;
  LJ_F_MIN[1] = LJ_SIGMA[1];
  LJ_F_MAX[1] = TWO_POW_ONE_SIX*LJ_SIGMA[1];
  B_MIN[1] = 0.32;
  B_MAX[1] = 0.42;
  LJ_RCSQ[1] = pow(3 * LJ_SIGMA[1],2);
  V_LJ_RC[1] = LJ_EPS[1] * (pow(ONE_THIRD, 12)-pow(ONE_THIRD,6));

  //type 2: carbon-hydrogen
  RC_SQ[2] = 3.24;
  REBO_F_MIN[2] = REBO_W_MIN[2] = 1.3;
  REBO_F_MAX[2] = 1.8;
  REBO_W_MAX[2] = 1.6;
  Q[2] = 0.340775728;
  A[2] = 149.94098723;
  ALPHA[2] = 4.10254983;
  B1[2] = 32.3551866587;
  B2[2] = 0.0;
  B3[2] = 0.0;
  BETA1[2] = 1.43445805925;
  BETA2[2] = 0.0;
  BETA3[2] = 0.0;
  LJ_SIGMA[2] = 3.025;
  LJ_EPS[2] = 4 * 0.002064;
  LJ_F_MIN[2] = LJ_SIGMA[2];
  LJ_F_MAX[2] = TWO_POW_ONE_SIX*LJ_SIGMA[2];
  B_MIN[2] = 0.75;
  B_MAX[2] = 0.90;
  LJ_RCSQ[2] = pow(3 * LJ_SIGMA[2],2);
  V_LJ_RC[2] = LJ_EPS[2] * (pow(ONE_THIRD, 12)-pow(ONE_THIRD,6));

  MASS[1] = 1.0079;
  MASS[6] = 12.011;
  MASS[18] = 39.948;
}
