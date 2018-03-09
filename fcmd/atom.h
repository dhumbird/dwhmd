//atom.h
//class listing for atom class

#ifndef ATOM_H
#define ATOM_H

#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>
#include <map>
#include "svector.h"
#include "chemtab.h"
#include "nbr.h"
#include "bicubic.h"
#include "string"

typedef vector<nbr>::iterator VNI;
extern map<short, double> MASS;

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
  //svector Fsave;
  double m, u, Nt;
  map<short, double> Nmap;
  vector<nbr> nlist;
  void * my_cell; void* p_clust;
  
  atom();
  void SetProps(){m=MASS[id]/9648.531;}
  void PosUpdate(double& dt, double& dtsq2, double& Lx, double& Ly){
    R+=dt*V+dtsq2*F/m;
    while (abs(R.x/Lx)>0.5) R.x += (R.x < 0) ? Lx : -Lx;
    while (abs(R.y/Ly)>0.5) R.y += (R.y < 0) ? Ly : -Ly;
  }
  void AddNbrs();
  VNI FindNbr(atom*);
  void EraseNbr(atom*);
  void ts_init() {
    F.clear(); 
    p_clust=NULL;
  }
  double Ek() {return 0.5*m*V.sqmag();}
  //void Sort_nlist();
  friend ostream& operator << (ostream&, atom&);
  friend istream& operator >> (istream&, atom&);
};

#endif
