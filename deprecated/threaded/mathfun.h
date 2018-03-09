//mathfun.h

#ifndef MATHFUN_H
#define MATHFUN_H

#include <cmath>
#include <cstdlib>

inline double SpExp(double arg) {return (arg<-50)? 0:exp(arg);}
inline double neg(double arg) {return -1*abs(arg);}
inline double rand01() {return(double)rand()/(double)RAND_MAX;}
#endif
