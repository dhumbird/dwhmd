//atom.h
//class listing for atom class

#ifndef ATOM_H
#define ATOM_H

#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>
#include "svector.h"
#include "chemtab.h"
#include "nbr.h"
#include "bicubic.h"
#include "ljnbr.h"

typedef vector<nbr>::iterator VNI;
typedef vector<ljnbr>::iterator VNJI;
class atom{

 public:

  //printed members
  svector R; //position
  svector V; //velocity
  int ix;
  short id;
  bool is_fixed;
  
  //non-printed members
  svector P, F;
  svector Fsave;
  double m, u, Nt, NH, NC;
  vector<nbr> nlist;
  vector<ljnbr> ljnlist;
  void * my_cell; void* p_clust;
  
  atom();
  void SetProps();
  void PosUpdate(double&, double&, double&, double&);
  void PostComp();
  void LJ_PreComp();
  VNI FindNbr(atom*);
  VNJI LJ_FindNbr(atom*);
  void EraseNbr(atom*);
  void LJ_EraseNbr(atom*);
  void ts_init() {F.clear(); Fsave.clear(); p_clust=NULL;}
  double Ek() {return 0.5*m*V.sqmag();}
  //void Sort_nlist();
  friend ostream& operator << (ostream&, atom&);
  friend istream& operator >> (istream&, atom&);
};

#endif
