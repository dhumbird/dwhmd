//c_vv.cpp
//member functions for the velocity verlet method

#include "config.h"

//************************************************************************
void config::FirstVV(){
  for (i=begin; i<end; i++){
    if (CheckSput(i)){
      if (!i->is_fixed){
	i->R+=dt*i->V+dtsq2*i->F/i->m;
	i->R.minimg(Lx, Ly, Lz);
	i->V+=dt2*i->F/i->m;
      }
    }
    else i--;
  }
}
//************************************************************************
void config::OFirstVV(){
  //FirstVV without CheckSput

  for (i=begin; i<end; i++){
    if (!i->is_fixed){
      i->R+=dt*i->V+dtsq2*i->F/i->m;
      i->R.minimg(Lx, Ly, Lz);
      i->V+=dt2*i->F/i->m;
    }
  }
}
//*******************************************************************
void config::SecondVV(){
  for (i=begin; i<end; i++)
    if (!i->is_fixed){
      i->V+=dt2*i->F/i->m;
      ek+= i->Ek();
    }
  T=ek/3/(N-NFixed)/KB;
  ek_av=ek/(N-NFixed);
}
