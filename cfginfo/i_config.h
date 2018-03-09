#ifndef I_CONFIG_H
#define I_CONFIG_H

#include <vector>
#include <set>
#include "i_atom.h"

typedef vector<atom>::iterator VAI;
typedef set<atom*>::iterator SAI;

#include <string>
#include <fstream>
#include <cstdio>
#include <cmath>
#include "consts.h"
#include "mathfun.h"
#include "strfun.h"
#include "sysfun.h"
#include <sstream>

//************************************************************************
class config{

 public:
  vector<atom>* atom_list;
  VAI begin, end;
  double t, Lx, Ly, Lz, T;
  int N, NFixed, Nmax;
  string name;
  double the_top;
  
  //******inline functions**********
  config(){init();}
  config(string scfg){init(); Load(scfg); ReNeighbor();}
  ~config(){delete atom_list;}
  void SetProps(){for (VAI i=begin; i<end; i++) i->SetProps();}
  double MinZ(){
    double z=1000; for (VAI w=begin; w<end; w++) z=mfmin(w->R.z,z); return z;}
  double MaxZ(){
    double z=-1000; for (VAI w=begin; w<end; w++) z=mfmax(w->R.z,z); return z;}

  void resize(int);
  VAI atomix(int);
  void init();
  void Load(string);
  void ReNeighbor();
};


#endif



