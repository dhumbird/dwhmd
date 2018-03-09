#include "ljnbr.h"

extern map<short,double> LJ_F_MIN, LJ_F_MAX;

//******************************************************************
void ljnbr::PreComp(){
  C = 1-C;
  r_true = sqrt(r_true);
  Rhat/=r_true;
  double ddmin, ddmax;
  ddmin=LJ_F_MIN[type]; ddmax=LJ_F_MAX[type];
  PolySwitch((r_true-ddmin)/(ddmax-ddmin), &f_LJ, &fprime_LJ);
  fprime_LJ/=(ddmax-ddmin);
}
