//atom.h
//class listing for atom class

#ifndef RATOM_H
#define RATOM_H

#include <cmath>
#include <vector>
#include <set>
#include <iostream>
#include "svector.h"
#include "rgb.h"

class atom{

 public:

  svector R; //position
  svector V; //velocity
  int ix;
  short id;
  bool is_fixed;
  
  //non-printed members
  rgb color;
  float radius; float alpha;
  double m;
  atom();
  set<atom*> nlist;
  double Ek() {return 0.5*m*V.sqmag();}
  void SetProps(float);
  friend ostream& operator << (ostream&, atom&);
  friend istream& operator >> (istream&, atom&);
};

#endif
