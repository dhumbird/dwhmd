#include "mathfun.h"

void PolySwitch(double t, double* S, double* Sprime){
  if (t < 0){
    *S=1;
    *Sprime=0;
  }
  else if (t < 1){
    *S = 1; *Sprime=0;
    double tt=t*t; //t^2
    *Sprime -= 30*tt;
    tt*=t; //t^3
    *S -= 10*tt;
    *Sprime += 60*tt;
    tt*=t; //t^4
    *S += 15*tt;
    *Sprime -=30*tt;
    *S-=6*tt*t;
  }
  else{
    *S=*Sprime=0;
  }
}
