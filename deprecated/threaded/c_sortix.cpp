#include "config.h"

class ix_less : public binary_function<particle,particle,bool>{
public:
  bool operator()(particle x, particle y){return x.ix < y.ix;}
};

class ix_eq : public binary_function<particle,int,bool>{
public:
  bool operator()(particle x, int index){return x.ix==index;}
};

//*********************************************************************
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




