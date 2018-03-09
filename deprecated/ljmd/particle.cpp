//particle.cpp
//function listing for particle class

#include "particle.h"
#include <math.h>

//***********constructors, destructor*******************
particle::particle(){
  R.clear();
  F.clear();
  V.clear();
  m=0.0;
  nlist=new short[1];
}


particle::~particle(){
  delete nlist;
}

void particle::listinit(int N){
  int i;
  delete nlist;
  nlist=new short[N];
  for (i=0; i<N; i++) nlist[i]=0;
}

ostream& operator << (ostream& out, particle& p){
  out<<p.m<<" "<<p.R<<" "<<p.V<<" "<<p.F;
  return out;
}

istream& operator >> (istream& in, particle& p){
  in>>p.m>>p.R>>p.V>>p.F;
  return in;
}





