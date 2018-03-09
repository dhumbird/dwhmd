#ifndef NBR_H
#define NBR_H

#include "mathfun.h"
#include "svector.h"
#include "consts.h"
#include "chemtab.h"
#include <string>

class nbr{
 public:
  nbr(){};
  void* a1;
  void* a2;
  double Fv_kl;
  double r, f, fprime;
  svector Rhat;
  void Invert();
  void PreComp();
  short type; 
  double bo;
};


#endif
