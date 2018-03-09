#include "config.h"
//*******************************************************************
void config::init(){
  ChemDef();
  atom=new vector<particle>;
  name="temp";
  t=u=ek=0; //doubles
  mat=0;
  Setdt(0.001);
  bhb_dt=0.01;
  trips_list.clear();
  Lcx=Lcy=Lcz=Lz_full=the_top=0; //doubles
  Nx=Ny=Nz=Nc=N=NFixed=Nmax=0; //ints
  Lx=Ly=Lz=T=0; //doubles
  cells.clear();
  sput_list.clear();
  inert=-1;
}
//**********************************************************************
void config::reset(){
  int index=0;
  for (i=begin; i<end; i++){
    i->ix=index++;
    i->O=i->R;
    i->xpass=i->ypass=0;
  }
  Nmax=N;
  t=0;
}



