#ifndef NBR_H
#define NBR_H

#include "mathfun.h"
#include "svector.h"
#include "consts.h"
#include "chemtab.h"

class nbr{
 public:
  nbr(){};
  void* a1;
  void* a2;
  double r, f, fprime, w, wprime;
  svector Rhat;
  void Invert();
  void PreComp();
  short type; //0 if C-C, 1 if H-H, 2 if C-H
};


#endif
