//atom.cpp

#include "ratom.h"

//***constructor**********************************************
atom::atom(){
  m=0;
  R.clear();
  V.clear();
  id=0; ix=0;
  is_fixed=0;
  alpha=0; radius=0;
}
//*************************************************************
void atom::SetProps(float radscale){
  switch (id){
  case 14:
    if (!is_fixed) color=rgb(0, 0, 1);
    else color=rgb(0, 0, 0.297);
    radius=1.11;
    m=28;
    break;
  case 9:
    color=rgb(1, 1, 1);
    radius=0.71;
    m=19;
    break;
  case 18:
    radius=0.97;
    color=rgb(0, 1, 0);
    m=40;
    break;
  case 6:
    radius=0.77;
    color=rgb(1, 0.5, 0);
    m=12;
    break;
  case 1:
    radius=0.79;
    color=rgb(1, 1, 1);
    m=1;
    break;
  case 17:
    radius=0.99;
    color=rgb(0, 1, 0);
    m=35;
    break;
  }
  radius*=radscale;
}
//**************************************************************
ostream& operator << (ostream& out, atom& p){
  out<<p.ix<<" "<<p.id<<" "<<p.R<<" "<<p.V<<" "<<p.is_fixed;
  return out;
}
//**************************************************************
istream& operator >> (istream& in, atom& p){
  in>>p.ix>>p.id>>p.R>>p.V>>p.is_fixed;
  return in;
}






