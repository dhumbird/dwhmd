//particle.cpp

#include "particle.h"

extern map<short, cheminfo> chemtab;

//***constructor**********************************************
particle::particle(){
  F.clear();
  nlist.clear();
  rc=m=0;
  three_body=0;
  R.clear();
  V.clear();
  P.clear();
  O.clear();
  id=0; ix=0;
  is_fixed=0;
  my_cell=p_clust=NULL;
  xpass=ypass=0;
}
//**************************************************************
void particle::SetProps(){
  m=chemtab[id].mass/9648.531;
  rc=chemtab[id].cutoff;
  three_body=chemtab[id].three_body;
}
//******************************************************************
void particle::PosUpdate(double& dt, double& dtsq2, double& Lx, double& Ly){
  R+=dt*V+dtsq2*F/m;
  if (abs(R.x/Lx)>0.5)
    if (R.x < 0){
      R.x+=Lx;
      xpass--;
    }
    else{
      R.x-=Lx;
      xpass++;
    }
  
  if (abs(R.y/Ly)>0.5)
    if (R.y < 0){
      R.y+=Ly;
      ypass--;
    }
    else{
      R.y-=Ly;
      ypass++;
    }
}
//**************************************************************
ostream& operator << (ostream& out, particle& p){
  out<<p.ix<<" "<<p.id<<" "<<p.R<<" "<<p.V<<" "<<p.O<<" "<<p.is_fixed<<" "<<
    p.xpass<<" "<<p.ypass;
  return out;
}
//**************************************************************
istream& operator >> (istream& in, particle& p){
  in>>p.ix>>p.id>>p.R>>p.V>>p.O>>p.is_fixed>>p.xpass>>p.ypass;
  return in;
}






