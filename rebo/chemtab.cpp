#include "chemtab.h"

map<short, double> MASS;
map<short,double> RC_SQ, REBO_F_MIN, REBO_F_MAX, Q, A, 
ALPHA, B1, B2, B3, BETA1, REBO_W_MIN, REBO_W_MAX; 
map<short,double>  BETA2, BETA3;

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

  MASS[1]=M_1;
  MASS[6]=M_6;
  MASS[18]=M_18;
}
