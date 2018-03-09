#include "config.h"

class Ek_less : public binary_function<particle, particle, bool>{
public:
  bool operator()(particle x, particle y){return x.V.sqmag() < y.V.sqmag();}
};

class ix_less : public binary_function<particle,particle,bool>{
public:
  bool operator()(particle x, particle y){return x.ix < y.ix;}
};

class ix_eq : public binary_function<particle,int,bool>{
public:
  bool operator()(particle x, int index){return x.ix==index;}
};

class Rz_less : public binary_function<particle, particle, bool>{
public:
  bool operator()(particle x, particle y){return  x.R.z < y.R.z;}
};
//**********************************************************************
void config::Ek_sort_asc(){
  sort(atom->begin(), atom->end(), Ek_less());
  begin=atom->begin(); end=atom->end();
}

void config::Ek_sort_desc(){
  Ek_sort_asc();
  reverse(atom->begin(), atom->end());
  begin=atom->begin(); end=atom->end();
}
//**********************************************************************
void config::ix_sort_asc(){
  sort(atom->begin(), atom->end(), ix_less());
  begin=atom->begin(); end=atom->end();
}

void config::ix_sort_desc(){
  ix_sort_asc();
  reverse(atom->begin(), atom->end());
  begin=atom->begin(); end=atom->end();
}
//********************************************************************
void config::Rz_sort_asc(){
  sort(atom->begin(), atom->end(), Rz_less());
  begin=atom->begin(); end=atom->end();
}

void config::Rz_sort_desc(){
  Rz_sort_asc();
  reverse(atom->begin(), atom->end());
  begin=atom->begin(); end=atom->end();
}





