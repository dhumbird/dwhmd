//c_atom.cpp
//config's member functions for manipulating individual atoms
#include "config.h"

VPI config::append(short id){
  resize(++N);
  VPI newatom=end-1;
  newatom->id=id;
  newatom->ix=Nmax++;
  newatom->SetProps();
  Partition(); ReNeighbor();
  return newatom;
}
//*************************************************************************
bool config::erase(VPI D){
  if (D==end) return 0;
  else{
    N--;
    atom->erase(D);
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
void config::resize(int s){
  N=s;
  atom->resize(s);
  begin=atom->begin();
  end=atom->end();
}

