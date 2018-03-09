#include "i_chemtab.h"

map<short,double> RC_SQ, F_MIN, F_MAX, MASS;

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
  RC_SQ[12] = 3.24;
  F_MIN[12] = 1.7;
  F_MAX[12] = 2.0;
  
  //silicon-silicon
  RC_SQ[28] = 7.29;
  F_MIN[28] = 2.7;
  F_MAX[28] = 3.0;
  
  //silicon-carbon
  F_MIN[20] = 2.204541;
  F_MAX[20] = 2.509980;
  RC_SQ[20] = pow(F_MIN[20],2);
  
  //silicon-fluorine
  F_MIN[23] = 1.83922;
  F_MAX[23] = 2.13922;
  RC_SQ[23] = pow(F_MIN[23],2);
  //fluorine-fluorine FROM TANAKA 2000
  RC_SQ[18] = 2.89;
  F_MIN[18] = 1.7;
  F_MAX[18] = 2.0;

  //carbon-fluorine FROM TANAKA 2000
  RC_SQ[15] = 2.89;
  F_MIN[15] = 1.7;
  F_MAX[15] = 2.0;
  //silicon-chlorine
  F_MIN[31] = 2.320;
  F_MAX[31] = 2.620;
  RC_SQ[31] = pow(F_MIN[31],2);
  //chlorine-chlorine
  F_MIN[34] = 2.284;
  F_MAX[34] = 2.584;
  RC_SQ[34] = pow(F_MIN[34],2);

}

