#include "config.h"

class Rz_less : public binary_function<particle, particle, bool>{
public:
  bool operator()(particle x, particle y){return  x.R.z < y.R.z;}
};

//  class Rx_less : public binary_function<particle, particle, bool>{
//  public:
//    bool operator()(particle x, particle y){return  x.R.x < y.R.x;}
//  };

//  class Ry_less : public binary_function<particle ,particle, bool>{
//  public:
//    bool operator()(particle x, particle y){return  x.R.y < y.R.y;}
//  };

//**********************************************************************
void config::Rz_sort_asc(){
  sort(atom->begin(), atom->end(), Rz_less());
  begin=atom->begin(); end=atom->end();
}

void config::Rz_sort_desc(){
  Rz_sort_asc();
  reverse(atom->begin(), atom->end());
  begin=atom->begin(); end=atom->end();
}
//**********************************************************************
//  void config::Rx_sort_asc(){
//    sort(atom->begin(), atom->end(), Rx_less());
//    begin=atom->begin(); end=atom->end();
//  }

//  void config::Rx_sort_desc(){
//    Rx_sort_asc();
//    reverse(atom->begin(), atom->end());
//    begin=atom->begin(); end=atom->end();
//  }
//  //**********************************************************************
//  void config::Ry_sort_asc(){
//    sort(atom->begin(), atom->end(), Ry_less());
//    begin=atom->begin(); end=atom->end();
//  }

//  void config::Ry_sort_desc(){
//    Ry_sort_asc();
//    reverse(atom->begin(), atom->end());
//    begin=atom->begin(); end=atom->end();
//  }
