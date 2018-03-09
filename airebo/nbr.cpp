#include "nbr.h"

extern map<short,double> REBO_F_MIN, REBO_F_MAX, REBO_W_MIN, REBO_W_MAX;

double dmin, dmax;
void nbr::Invert(){
  void* dum=a1;
  a1=a2;
  a2=dum;
  Rhat*=-1;
}

//******************************************************************
void nbr::PreComp(){
  dmin=REBO_F_MIN[type]; dmax=REBO_F_MAX[type];
  PolySwitch((r-dmin)/(dmax-dmin), &f, &fprime);
  fprime/=(dmax-dmin);
  
  if (type == 2){
    dmin=REBO_W_MIN[type]; dmax=REBO_W_MAX[type];
    PolySwitch((r-dmin)/(dmax-dmin), &w, &wprime);
    wprime/=(dmax-dmin);
  }
  else{
    w=f;
    wprime=fprime;
  }
}
