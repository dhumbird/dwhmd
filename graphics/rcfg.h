#ifndef RCFG_H
#define RCFG_H

#include <vector>
#include <set>
#include "ratom.h"
#include <iostream>
#include <sstream>
#include <map>
#include "miscfun.h"
#include <algorithm>
using namespace std;
typedef vector<atom>::iterator VAI;
typedef set<atom*>::iterator SAI;

//************************************************************************
class rcfg{
  
 public:
  vector<atom>* atom_list;
  VAI begin, end, i, j, k, l;
  double t, Lx, Ly, Lz;
  int N, Nmax;
  string name;
  SAI sai_i, sai_j;
  double the_top, the_bottom;
  //******inline functions**********
  rcfg(){init();}
  rcfg(string scfg, float rad){init(); Load(scfg, rad); ReNeighbor();}
  double MinZ(){
    double z=1000; for (VAI w=begin; w<end; w++) z=w->R.z<?z; return z;}
  double MaxZ(){
    double z=-1000; for (VAI w=begin; w<end; w++) z=w->R.z>?z; return z;}
  //******in rcfg.cpp************
  void ReNeighbor();
  void resize(int);
  VAI atomix(int);
  void init();
  void Load(string,float);
  void ix_sort_asc();
  void ix_sort_desc();
};


#endif



