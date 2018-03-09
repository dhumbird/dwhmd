//mathfun.h

#ifndef MATHFUN_H
#define MATHFUN_H

#include <cmath>
#include <cstdlib>
#include <ctime>
inline double SpExp(double arg) {return (arg<-50)? 0:exp(arg);}
inline double neg(double arg) {return -1*abs(arg);}
inline double rand01() {return(double)rand()/(double)RAND_MAX;}
inline int TimeSeed(){return (unsigned)clock()+(unsigned)time(NULL);}
#endif
