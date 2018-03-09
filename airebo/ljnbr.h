#ifndef LJNBR_H
#define LJNBR_H

#include "nbr.h"

class ljnbr : public nbr{
 public:
  ljnbr(){
    f_LJ=fprime_LJ=C=0;
    n1=n2=n3=NULL;
  };
/*   void* a1; */
/*   void* a2; */
/*   double r, f, fprime, w, wprime; */
  double r_true;
  double f_LJ, fprime_LJ, C;
  nbr *n1, *n2, *n3;
/*   svector Rhat; */
  void PreComp(); 
/*   short type; //0 if C-C, 1 if H-H, 2 if C-H */
};


#endif
