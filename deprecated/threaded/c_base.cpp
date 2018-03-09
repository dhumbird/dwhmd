#include "config.h"

//*******************************************************************
config::config(){
  init();
}

config::config(string scfg){
  init();
  Load(scfg);
  Partition();
}
//*******************************************************************
void config::init(){
  chem_is_def=0;
  atom=new vector<particle>;
  name="temp";
  t=0; u=0; ek=0;
  mat=0;
  Setdt(0.001);
  bhb_dt=0.01;
  a=0; d1=0; Lcx=0; Lcy=0; Lcz=0; sLz=0; sLz2=0;
  Nx=0; Ny=0; Nz=0; Nc=0; N=0; NFixed=0; Nsput=0;
  Lx=0; Ly=0; Lz=0; T=0; Lz2=0; ek_av=0; u_av=0;
  cells.clear();
}
//********************************************************************
void config::resize(int s){
  N=s;
  atom->resize(s);
  begin=atom->begin();
  end=atom->end();
}




