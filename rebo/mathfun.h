//mathfun.h

#ifndef MATHFUN_H
#define MATHFUN_H

#include <cmath>
#include <cstdlib>
#include <ctime>
#include "consts.h"

inline double SpExp(double arg) {return (arg<-50)? 0:exp(arg);}
inline double neg(double arg) {return -1*fabs(arg);}
inline double rand01() {return(double)rand()/(double)RAND_MAX;}
inline int TimeSeed(){return (unsigned)clock()+(unsigned)time(NULL);}
void Gpoly(const double, int, double*, double*, double*, double*);
inline int signof(double d){return d < 0 ? -1:1;}
inline double mfmin(double a, double b){
  if (a<b) return a; else return b;
}
inline double mfmax(double a, double b){
  if (a>b) return a; else return b;
}
void PolySwitch(double, double*, double*);
#endif
