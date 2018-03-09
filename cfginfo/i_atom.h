//atom.h
//class listing for atom class

#ifndef I_ATOM_H
#define I_ATOM_H

#include <cmath>
#include <set>
#include <iostream>
#include "svector.h"
#include "i_chemtab.h"
using namespace std;

class atom{
 public:
  svector R; //position
  svector V; //velocity
  int ix;
  short id;
  bool is_fixed;
  double m;
  set<atom*> nlist;
    
  atom();
  void SetProps();
  double Ek() {return 0.5*m*V.sqmag();}
  friend ostream& operator << (ostream&, atom&);
  friend istream& operator >> (istream&, atom&);
};

#endif
