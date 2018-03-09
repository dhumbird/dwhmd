#include "config.h"

class Ek_less : public binary_function<particle, particle, bool>{
public:
  bool operator()(particle x, particle y){return x.V.sqmag() < y.V.sqmag();}
};

class F_less : public binary_function<particle, particle, bool>{
public:
  bool operator()(particle x, particle y){return x.F.sqmag() < y.F.sqmag();}
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
void config::F_sort_asc(){
  sort(atom->begin(), atom->end(), F_less());
  begin=atom->begin(); end=atom->end();
}

void config::F_sort_desc(){
  F_sort_asc();
  reverse(atom->begin(), atom->end());
  begin=atom->begin(); end=atom->end();
}
