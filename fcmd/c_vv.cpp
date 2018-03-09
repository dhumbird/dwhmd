//c_vv.cpp
//member functions for the velocity verlet method

#include "config.h"

//************************************************************************
void config::FirstVV(){
  for (VAI i=begin; i<end; i++){
    if (!i->is_fixed){
      i->PosUpdate(dt, dtsq2, Lx, Ly);
      i->V+=dt2*i->F/i->m;
      if ((this_cell=WhichCell(i->R)) != i->my_cell){
        ((subcell*) i->my_cell)->erase(&(*i));
        ((subcell*) this_cell)->insert(&(*i));
      }
    }
  }
}

//*******************************************************************
void config::SecondVV(){
  for (VAI i=begin; i<end; i++)
    if (!i->is_fixed){
      i->V+=dt2*i->F/i->m;
      ek+= i->Ek();
    }
  T=ek/3/(float)(N-NFixed)/KB;
}


