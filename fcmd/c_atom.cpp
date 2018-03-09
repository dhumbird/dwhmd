//c_atom.cpp
//config's member functions for manipulating individual atoms
#include "config.h"

//*************************************************************************
bool config::erase(VAI D){
  if (D==end) return 0;
  else{
    N--;
    if (D->ix==Nmax-1) Nmax--;
    atom_list->erase(D);
    begin=atom_list->begin();
    end=atom_list->end();
    Partition(); ReNeighbor();
    return 1;
  }
}
//*****************************************************************************
VAI config::atomix(int index){
  VAI p=begin;
  while (p<end && p->ix!=index) p++;
  return p;
}
//**************************************************************************
void config::resize(int s){
  N=s;
  atom_list->resize(s);
  begin=atom_list->begin();
  end=atom_list->end();
}
//****************************************************************************
VAI config::append(short id){
  resize(++N);
  VAI newatom=end-1;
  newatom->id=id;
  newatom->ix=Nmax++;
  newatom->SetProps();
  Partition(); ReNeighbor();
  return newatom;
}


