//particle.h
//class listing for particle class

#ifndef PARTICLE_H
#define PARTICLE_H

#include <cmath>
#include <set>
#include <iostream>
#include "svector.h"
#include "chemtab.h"

class particle{

 public:
  svector R; //position
  svector P; //multi-purpose extra svector
  svector V; //velocity
  int ix; //index in list
  short id; //atomic number; represents species
  double m;
  double u;
  bool is_fixed;
  svector F; //force
  set<particle*> nlist;
  double rc; //squared cutoff, sq. ang
  bool three_body;
  particle();
  void SetProps();
  inline void ts_init(){u=0; nlist.clear(); F.clear();}
  void Copy(particle&);
  inline double Ek() {return 0.5*m*V.sqmag();}
  friend ostream& operator << (ostream&, particle&);
  friend istream& operator >> (istream&, particle&);
};

#endif
