#include "config.h"

//**********************************************************************
void config::reset(){
  int index=0;
  for (i=begin; i<end; i++){
    i->ix=index++;
  }
  Nmax=N;
  t=0;
}
//*********************************************************************
void config::init(){
  ChemDef();
  atom_list=new vector<atom>;
  name="temp";
  t=u=ek=0; //doubles
  Setdt(0.0001);
  inert=-1;
  N=NFixed=Nmax=0; //ints
  Lx=Ly=Lz=T=0; //doubles
}


