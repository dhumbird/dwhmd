#include "nbr.h"

extern map<short,double> F_MIN, F_MAX;

double dmin, dmax;
void nbr::Invert(){
  void* dum=a1;
  a1=a2;
  a2=dum;
  Rhat*=-1;
}

//******************************************************************
void nbr::PreComp(){
  dmin=F_MIN[type]; dmax=F_MAX[type];
  PolySwitch((r-dmin)/(dmax-dmin), &f, &fprime);
  fprime/=(dmax-dmin);
  bo=0;
}
