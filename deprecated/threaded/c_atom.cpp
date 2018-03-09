//c_atom.cpp
//config's member functions for manipulating individual atoms
#include "config.h"
//************************************************************************
void config::CopyAtom(VPI C){
  resize(++N);
  (end-1)->Copy(*C);
  RePartition();
}
//*************************************************************************
bool config::erase(VPI D){
  if (D==end) return 0;
  else{
    atom->erase(D);
    N--;
    begin=atom->begin();
    end=atom->end();
    return 1;
  }
}
//*****************************************************************************
VPI config::atomix(int index){
  p=begin;
  while (p<end && p->ix!=index) p++;
  return p;
}
//**************************************************************************
