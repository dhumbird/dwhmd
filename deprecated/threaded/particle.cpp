//particle.cpp

#include "particle.h"

extern map<short, cheminfo> chemtab;

//***constructor**********************************************
particle::particle(){
  F.clear();
  nlist.clear();
  rc=0;
  u=0;
  three_body=0;
  R.clear();
  V.clear();
  P.clear();
  id=0; ix=0;
  is_fixed=0;
  m=0;
}
//**************************************************************
void particle::SetProps(){
  m=chemtab[id].mass/9648.531;
  rc=chemtab[id].cutoff;
  three_body=chemtab[id].three_body;
}
//**************************************************************
void particle::Copy(particle& p){
  R=p.R;
  P=p.P;
  V=p.V;
  F=p.F;
  ix=0;
  id=p.id;
  SetProps();
  is_fixed=p.is_fixed;
}
//**************************************************************
ostream& operator << (ostream& out, particle& p){
  out<<p.ix<<" "<<p.id<<" "<<p.R<<" "<<p.V<<" "<<p.is_fixed;
  return out;
}

//**************************************************************
istream& operator >> (istream& in, particle& p){
  in>>p.ix>>p.id>>p.R>>p.V>>p.is_fixed;
  return in;
}






