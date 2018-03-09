//atom.cpp

#include "i_atom.h"

extern map<short, double> MASS;

//***constructor**********************************************
atom::atom(){
  nlist.clear();
  m=0;
  R.clear();
  V.clear();
  id=0; ix=0;
  is_fixed=0;
}
//**************************************************************
void atom::SetProps(){
  m=MASS[id]/9648.531;
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






